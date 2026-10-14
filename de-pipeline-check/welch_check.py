from scipy import stats
from statsmodels.stats.multitest import multipletests
import numpy as np

def welch_ttest_de(data, protein_ids, control_idx, treatment_idx, log2fc_threshold=1.0, p_threshold=0.05):
    results = []
    p_values = []

    for i in range(data.shape[0]):
        control = data[i, control_idx]
        treatment = data[i, treatment_idx]

        # Remove NaN
        control = control[~np.isnan(control)]
        treatment = treatment[~np.isnan(treatment)]

        if len(control) < 2 or len(treatment) < 2:
            continue

        log2fc = np.mean(treatment) - np.mean(control)  # Already log2
        _, pval = stats.ttest_ind(control, treatment, equal_var=False)

        results.append({
            'protein': protein_ids[i],
            'log2FC': log2fc,
            'pvalue': pval,
        })
        p_values.append(pval)

    # BH correction (CRITICAL — never skip this!)
    _, adj_pvals, _, _ = multipletests(p_values, method='fdr_bh')
    for i, r in enumerate(results):
        r['adj_pvalue'] = adj_pvals[i]
        r['significant'] = abs(r['log2FC']) >= log2fc_threshold and r['adj_pvalue'] < p_threshold   # threshold on the log2 scale

    return results


rng=np.random.default_rng(0); data=rng.normal(20,0.3,(200,6)); data[:20,3:]+=2; data[5,0]=np.nan
protein_ids=[f'P{i}' for i in range(200)]
res=welch_ttest_de(data, protein_ids, [0,1,2], [3,4,5])
print('tested', len(res), 'significant', sum(r['significant'] for r in res))

