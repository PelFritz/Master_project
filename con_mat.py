import pandas as pd
import seaborn as sn
import matplotlib.pyplot as plt

# Model(96 filters) trained on promoters with snp only
con_mat1 = [[10, 33662],
            [1847, 606871]]
df_cm1 = pd.DataFrame(con_mat1, range(2), range(2))
sn.set(font_scale=1.4)
sn.heatmap(df_cm1, annot=True, annot_kws={'size': 16}, cmap='coolwarm', fmt='g')
plt.xlabel('Predicted')
plt.ylabel('True')
plt.title('Confusion matrix')
plt.show()

# Model(64 filters) trained on promoters with snp only
con_mat2 = [[10, 33662],
            [2114, 606604]]
df_cm2 = pd.DataFrame(con_mat2, range(2), range(2))
sn.set(font_scale=1.4)
sn.heatmap(df_cm2, annot=True, annot_kws={'size': 16}, cmap='coolwarm', fmt='g')
plt.xlabel('Predicted')
plt.ylabel('True')
plt.title('Confusion matrix')
plt.show()

# Model (64 filters) trained on promoters with snp only and tested on those with snp + indels
con_mat3 = [[72, 33574],
            [2553, 605225]]
df_cm3 = pd.DataFrame(con_mat3, range(2), range(2))
sn.set(font_scale=1.4)
sn.heatmap(df_cm3, annot=True, annot_kws={'size': 16}, cmap='coolwarm', fmt='g')
plt.xlabel('Predicted')
plt.ylabel('True')
plt.title('Confusion matrix')
plt.show()

# Model output for scenario B
con_mat4 = [[1121, 33525],
            [16935, 590843]]
df_cm4 = pd.DataFrame(con_mat4, range(2), range(2))
sn.set(font_scale=1.4)
sn.heatmap(df_cm4, annot=True, annot_kws={'size': 16}, cmap='coolwarm', fmt='g')
plt.xlabel('Predicted')
plt.ylabel('True')
plt.title('Confusion matrix')
plt.show()

# Model output for scenario C

con_mat5 = [[16122, 88944],
            [40096, 496262]]
df_cm5 = pd.DataFrame(con_mat5, range(2), range(2))
sn.set(font_scale=1.4)
sn.heatmap(df_cm5, annot=True, annot_kws={'size': 16}, cmap='coolwarm', fmt='g')
plt.xlabel('Predicted')
plt.ylabel('True')
plt.title('Confusion matrix')
plt.show()
