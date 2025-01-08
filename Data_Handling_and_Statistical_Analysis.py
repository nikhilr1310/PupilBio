# Import necessary libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import variation, fisher_exact
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, confusion_matrix

# Settings for reproducibility
import random
np.random.seed(42)
random.seed(42)

# Load the dataset
df = pd.read_csv('PupilBioTest_PMP_revA.csv')
print("Columns in the dataset:")
print(df.columns)

# Data handling
print(df.info())
print(df.describe())

'''a)Coverage Analysis:
    a. Calculate the median and coefficient of variation (CV) for single CpG coverage
       in each tissue '''

# Calculate total coverage for each row
status_columns = ['`000', '`001', '`010', '`011', '`100', '`101', '`110', '`111']
df['Coverage'] = df[status_columns].sum(axis=1)

# Calculate the median and coefficient of variation separately
median_coverage = df.groupby('Tissue')['Coverage'].median().reset_index(name='Median')
std_coverage = df.groupby('Tissue')['Coverage'].std().reset_index(name='Std')
mean_coverage = df.groupby('Tissue')['Coverage'].mean().reset_index(name='Mean')

# Merge statistics into one DataFrame
coverage_stats = pd.merge(median_coverage, mean_coverage, on='Tissue')
coverage_stats = pd.merge(coverage_stats, std_coverage, on='Tissue')

# Calculate coefficient of variation manually
coverage_stats['coefficient_of_variation'] = coverage_stats['Std'] / coverage_stats['Mean']*100

print("Coverage Statistics (Median, Mean, Std, and coefficient of variation):")
print(coverage_stats)

#b. Generate plots summarizing the coverage statistics
print("Columns in coverage_stats DataFrame:")
print(coverage_stats.columns)

# Boxplot for Single CpG Coverage Distribution by Tissue
plt.figure(figsize=(10, 6))
sns.boxplot(x='Tissue', y='Coverage', data=df, hue='Tissue', dodge=False, palette='Set2')
plt.title('Single CpG Coverage Distribution by Tissue')
plt.xlabel('Tissue')
plt.ylabel('Coverage')
plt.legend([], [], frameon=False)  # Remove redundant legend
# plt.show()

# Bar plot for CV comparison
plt.figure(figsize=(8, 5))
sns.barplot(x='Tissue', y='coefficient_of_variation', data=coverage_stats, hue='Tissue', dodge=False, palette='Set2')
plt.title('Coefficient of Variation (CV) for Single CpG Coverage')
plt.xlabel('Tissue')
plt.ylabel('Coefficient of Variation')
plt.legend([], [], frameon=False)  # Remove redundant legend
# plt.show()

'''Biomarker Identification: 
    a. Identify PMPs with high specificity for tissue differentiation, minimizing false 
       positives for Tissue #1 while allowing some false negatives. Use statistical or 
       machine learning approaches to assign confidence (e.g., p-values) to each PMP'''

# Step 1: Calculate PMP occurrences
status_columns = ['`000', '`001', '`010', '`011', '`100', '`101', '`110', '`111']
df['PMP'] = df['CpG_Coordinates'] + "_" + df['strand']  # Create unique PMP identifier

pmp_counts = df.groupby(['PMP', 'Tissue'])[status_columns].sum().reset_index()
pmp_counts['Tissue_1_Count'] = pmp_counts['`000'] + pmp_counts['`001']  # Adjust based on actual count distribution
pmp_counts['Tissue_2_Count'] = pmp_counts['`110'] + pmp_counts['`111']  # Adjust based on actual count distribution
pmp_counts['Total_Count'] = pmp_counts['Tissue_1_Count'] + pmp_counts['Tissue_2_Count']

# Step 2: Calculate specificity and sensitivity
pmp_counts['Specificity'] = pmp_counts['Tissue_2_Count'] / pmp_counts['Total_Count']
pmp_counts['Sensitivity'] = pmp_counts['Tissue_2_Count'] / (
    pmp_counts['Tissue_2_Count'] + pmp_counts['Tissue_1_Count'])

# Step 3: Assign p-values using Fisher's Exact Test
def calculate_p_value(row):
    contingency_table = [
        [row['Tissue_1_Count'], row['Tissue_2_Count']],
        [pmp_counts['Tissue_1_Count'].sum(), pmp_counts['Tissue_2_Count'].sum()]
    ]
    _, p_value = fisher_exact(contingency_table, alternative='greater')
    return p_value

pmp_counts['p_value'] = pmp_counts.apply(calculate_p_value, axis=1)

# Filter PMPs with high specificity and significant p-values
filtered_pmps = pmp_counts[
    (pmp_counts['Specificity'] > 0.9) & (pmp_counts['p_value'] < 0.05)
]

# Display top PMPs
print("High Specificity PMPs:")
print(filtered_pmps.head())

# Step 4: Machine Learning for Classification
X = df[status_columns]  # Features: methylation status patterns
y = df['Tissue']        # Target: Tissue type

# Train-Test Split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, stratify=y, random_state=42)

# Random Forest Classifier
rf_model = RandomForestClassifier(random_state=42)
rf_model.fit(X_train, y_train)

# Predictions
y_pred = rf_model.predict(X_test)

# Classification Report
print("Classification Report:")
print(classification_report(y_test, y_pred))

# Feature Importance
feature_importance = pd.DataFrame({
    'Feature': X.columns,
    'Importance': rf_model.feature_importances_
}).sort_values(by='Importance', ascending=False)

# Display Top Features
print("Top Features (PMP Importance):")
print(feature_importance.head())

# Step 5: Visualization
plt.figure(figsize=(10, 6))
sns.barplot(data=feature_importance.head(10), x='Importance', y='Feature', palette='coolwarm',hue=None,legend=False)
plt.title("Top 10 PMP Features by Importance")
plt.xlabel("Feature Importance")
plt.ylabel("PMP")
plt.tight_layout()
plt.show()


# Step 1: Aggregate Reads for Each PMP and Tissue
status_columns = ['`000', '`001', '`010', '`011', '`100', '`101', '`110', '`111']

# Sum methylation patterns for each PMP and Tissue
pmp_counts = df.groupby(['PMP', 'Tissue'])[status_columns].sum().reset_index()

# Step 2: Calculate Total Reads and Variant Read Fraction (VRF)
# Convert counts to integer for arithmetic operations
pmp_counts[status_columns] = pmp_counts[status_columns].astype(int)

# Sum to get total reads
pmp_counts['Total_Reads'] = pmp_counts[status_columns].sum(axis=1)

# Sum of variant methylation patterns
pmp_counts['Variant_Reads'] = pmp_counts[['`001', '`010', '`011', '`100', '`101', '`110', '`111']].sum(axis=1)

# Mean Variant Read Fraction (VRF) calculation
pmp_counts['VRF'] = pmp_counts['Variant_Reads'] / pmp_counts['Total_Reads']

print("Mean Variant Read Fraction (VRF) for each PMP:")
print(pmp_counts[['PMP', 'Tissue', 'VRF']].head())