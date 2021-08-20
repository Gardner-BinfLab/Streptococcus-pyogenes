#### Jupyter notebooks
```41588_2019_417_MOESM3_ESM.csv``` is Supplementary Table 2 from [Davies et al., 2019](https://doi.org/10.1038/s41588-019-0417-8).
`classifiers.ipynb` is the correct file to use. Everything 
else are temporary WIP scripts. 


It includes AUCs (Area under the ROC curve) of the bitscores from eggNOG. The 
classification follows as 1 vs all scheme. For example: for the case `Superficial soft tissue infection`, this is 
labelled as 1, everything else is labelled as 0. Based on max AUC scores, we then pick the top 10 eggNOG models for 
each case:

 ```
 {'Superficial soft tissue infection': Index(['304Z7.faa.final_tree.fa', 'COG3649.faa.final_tree.fa',
        '33BC7.faa.final_tree.fa', '32R25.faa.final_tree.fa',
        'COG3646.faa.final_tree.fa', 'COG1312.faa.final_tree.fa',
        'COG1674.faa.final_tree.fa', '330BC.faa.final_tree.fa',
        '32UC6.faa.final_tree.fa', '33RTZ.faa.final_tree.fa'],
       dtype='object'),
 'iGAS': Index(['COG3649.faa.final_tree.f....
 ```


These highest scoring 10 models are then used to make random forest classifiers for classification of different cases. 
