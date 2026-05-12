import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split, StratifiedKFold, GridSearchCV
from sklearn.metrics import (accuracy_score, precision_score, recall_score, f1_score,
                             classification_report, confusion_matrix, precision_recall_curve,
                             average_precision_score, roc_auc_score)
from sklearn.preprocessing import LabelEncoder
from xgboost import XGBClassifier
import pickle
import os
import shap

shap.initjs()

def preprocess_data(file_path, target_column, test_size=0.2, random_state=42):
    data = pd.read_csv(file_path)
    X = data.drop(columns=[target_column])
    y = data[target_column]

    label_mapping = None
    if y.dtype == 'object' or y.dtype.name == 'category':
        le = LabelEncoder()
        y = le.fit_transform(y)
        label_mapping = dict(zip(le.classes_, le.transform(le.classes_)))

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=test_size, stratify=y, random_state=random_state
    )

    X_train = X_train.reset_index(drop=True)
    X_test = X_test.reset_index(drop=True)
    y_train = pd.Series(y_train, name=target_column).reset_index(drop=True)
    y_test = pd.Series(y_test, name=target_column).reset_index(drop=True)

    class_counts = dict(zip(*np.unique(y_train, return_counts=True)))
    print(f"\nClass distribution: {class_counts}")
    scale_pos_weight = class_counts[0] / class_counts[1] if 1 in class_counts else 1

    return X_train, X_test, y_train, y_test, label_mapping, scale_pos_weight

def calculate_metrics(y_true, y_pred, y_proba):
    return {
        "accuracy": float(accuracy_score(y_true, y_pred)),
        "precision": float(precision_score(y_true, y_pred)),
        "recall": float(recall_score(y_true, y_pred)),
        "f1": float(f1_score(y_true, y_pred)),
        "auc": float(roc_auc_score(y_true, y_proba)),
        "pr_auc": float(average_precision_score(y_true, y_proba)),
        "classification_report": classification_report(y_true, y_pred, output_dict=True),
        "y_true": list(map(int, y_true)),
        "y_proba": list(map(float, y_proba)),
        "y_pred": list(map(int, y_pred))
    }

def train_xgboost_fast(X_train, X_test, y_train, y_test, param_grid, scale_pos_weight, random_state=42):
    X_train_np = X_train.values
    y_train_np = y_train.values
    X_test_np = X_test.values
    y_test_np = y_test.values

    base_params = {
        'random_state': random_state,
        'eval_metric': 'logloss',
        'scale_pos_weight': scale_pos_weight
    }

    grid_search = GridSearchCV(
        estimator=XGBClassifier(**base_params),
        param_grid=param_grid,
        scoring={'auc': 'roc_auc', 
                 'precision': 'precision',
                 'recall': 'recall',
                 'f1': 'f1'},
        refit='auc',
        cv=StratifiedKFold(n_splits=10, shuffle=True, random_state=random_state),
        n_jobs=-1,
        verbose=1
    )
    grid_search.fit(X_train_np, y_train_np)

    best_params = grid_search.best_params_
    print(f"Best parameters: {best_params}")

    final_model = grid_search.best_estimator_
    final_model.fit(X_train_np, y_train_np)

    test_proba = final_model.predict_proba(X_test_np)[:, 1]
    test_pred = final_model.predict(X_test_np)
    test_metrics = calculate_metrics(y_test_np, test_pred, test_proba)

    with open("final_model.pkl", "wb") as f:
        pickle.dump(final_model, f)

    explainer = shap.TreeExplainer(final_model)
    shap_values = explainer.shap_values(X_train_np)

    cv_metrics = {
        metric: {
            'mean': grid_search.cv_results_[f'mean_test_{metric}'][grid_search.best_index_],
            'std': grid_search.cv_results_[f'std_test_{metric}'][grid_search.best_index_]
        } for metric in ['auc', 'precision', 'recall', 'f1']
    }

    return grid_search.cv_results_, best_params, test_metrics, shap_values, final_model, cv_metrics

class Visualizer:
    def __init__(self, output_path):
        self.output_path = output_path
        os.makedirs(output_path, exist_ok=True)

    def save_plot_data(self, data, filename):
        if isinstance(data, dict):
            df = pd.DataFrame([data])
        elif isinstance(data, pd.DataFrame):
            df = data
        else:
            df = pd.DataFrame(data)
        df.to_csv(f"{self.output_path}/{filename}_data.csv", index=False)

    def plot_combined_metrics(self, cv_metrics, test_metrics):
        metrics = ['auc', 'precision', 'recall', 'f1']
        df = pd.DataFrame([{
            'Metric': metric.upper(),
            'CV Mean': cv_metrics[metric]['mean'],
            'CV Std': cv_metrics[metric]['std'],
            'Test Value': test_metrics[metric]
        } for metric in metrics])
        
        self.save_plot_data(df, "combined_metrics")

        plt.figure(figsize=(10, 6))
        x = np.arange(len(metrics))
        width = 0.35
        plt.bar(x - width/2, df['CV Mean'], yerr=df['CV Std'],
                width=width, label='CV Mean ± Std', alpha=0.8)
        plt.bar(x + width/2, df['Test Value'], width=width, 
                label='Test Set', alpha=0.8)
        
        plt.xticks(x, df['Metric'])
        plt.ylabel('Score')
        plt.ylim(0, 1.1)
        plt.title('Performance Metrics Comparison')
        plt.legend()
        plt.tight_layout()
        plt.savefig(f"{self.output_path}/combined_metrics.pdf")
        plt.close()

    def plot_pr_curve(self, test_metrics):
        precision, recall, _ = precision_recall_curve(
            test_metrics['y_true'], test_metrics['y_proba']
        )
        test_ap = average_precision_score(test_metrics['y_true'], test_metrics['y_proba'])

        pr_data = pd.DataFrame({"recall": recall, "precision": precision})
        self.save_plot_data(pr_data, "pr_curve")

        plt.figure(figsize=(8, 6))
        plt.plot(recall, precision, color='red', linestyle='--',
                 label=f'Test Set (AP = {test_ap:.2f})')
        plt.xlabel('Recall')
        plt.ylabel('Precision')
        plt.title('Precision-Recall Curve')
        plt.legend()
        plt.tight_layout()
        plt.savefig(f"{self.output_path}/pr_curve.pdf")
        plt.close()

    def plot_confusion_matrix(self, test_metrics, class_names):
        cm = confusion_matrix(test_metrics['y_true'], test_metrics['y_pred'])
        cm_df = pd.DataFrame(cm, index=class_names, columns=class_names)
        self.save_plot_data(cm_df, "confusion_matrix")

        plt.figure(figsize=(6, 6))
        sns.heatmap(cm, annot=True, fmt='d', cmap='Blues',
                    xticklabels=class_names, yticklabels=class_names)
        plt.xlabel('Predicted')
        plt.ylabel('Actual')
        plt.title('Confusion Matrix')
        plt.savefig(f"{self.output_path}/confusion_matrix.pdf")
        plt.close()

    def plot_shap_summary(self, shap_values, X_train):
        plt.figure()
        shap.summary_plot(shap_values, X_train, feature_names=X_train.columns, show=False)
        plt.tight_layout()
        plt.savefig(f"{self.output_path}/shap_summary.pdf")
        plt.close()

    def plot_feature_importance(self, model, feature_names):
        importances = model.feature_importances_
        indices = np.argsort(importances)[::-1]
        sorted_features = [feature_names[i] for i in indices]
        sorted_importances = importances[indices]
        
        df = pd.DataFrame({'Feature': sorted_features, 'Importance': sorted_importances})
        self.save_plot_data(df, "feature_importance")
        
        plt.figure(figsize=(10, 6))
        plt.barh(range(len(sorted_features)), sorted_importances, align='center')
        plt.yticks(range(len(sorted_features)), sorted_features)
        plt.xlabel('Importance Score')
        plt.title('XGBoost Feature Importance')
        plt.gca().invert_yaxis()
        plt.tight_layout()
        plt.savefig(f"{self.output_path}/feature_importance.pdf")
        plt.close()

def main():
    parser = argparse.ArgumentParser(description='Optimized XGBoost Classifier')
    parser.add_argument('--data', required=True, help='Input CSV file path')
    parser.add_argument('--target', required=True, help='Target column name')
    parser.add_argument('--output_path', default='results', help='Output directory')
    args = parser.parse_args()

    viz = Visualizer(args.output_path)

    X_train, X_test, y_train, y_test, label_map, scale_pos_weight = preprocess_data(
        args.data, args.target
    )

    param_grid = {
        'n_estimators': [100, 200],
        'max_depth': [3, 5],
        'learning_rate': [0.05, 0.1],
        'subsample': [0.6, 0.8],
        'colsample_bytree': [0.6, 0.8],
        'reg_alpha': [0, 0.1],
        'reg_lambda': [0, 0.1]
    }

    cv_results, best_params, test_metrics, shap_values, final_model, cv_metrics = train_xgboost_fast(
        X_train, X_test, y_train.astype(int), y_test.astype(int),
        param_grid, scale_pos_weight
    )

    metadata = {
        'best_params': {k: (v.item() if isinstance(v, np.generic) else v) for k, v in best_params.items()},
        'test_metrics': {k: (v.item() if isinstance(v, np.generic) else v) if not isinstance(v, dict) else v
                         for k, v in test_metrics.items()},
        'class_distribution': {int(k): int(v) for k, v in zip(*np.unique(y_train, return_counts=True))},
        'scale_pos_weight': float(scale_pos_weight),
        'feature_names': list(X_train.columns)
    }
    pd.DataFrame([metadata]).to_csv(f'{args.output_path}/metadata.csv', index=False)

    class_names = list(label_map.keys()) if label_map else ['Negative', 'Positive']
    viz.plot_combined_metrics(cv_metrics, test_metrics)
    viz.plot_pr_curve(test_metrics)
    viz.plot_confusion_matrix(test_metrics, class_names)
    viz.plot_shap_summary(shap_values, X_train)
    viz.plot_feature_importance(final_model, X_train.columns)

if __name__ == '__main__':
    main()
