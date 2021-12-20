import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def _melt_table(
	abundance_table: pd.DataFrame,
	metadata: pd.DataFrame
) -> pd.DataFrame:
	# melt the table so that we can merge it with metadata later
	melt = abundance_table.unstack()
	melt_table = pd.DataFrame(melt)
	melt_table.reset_index(inplace=True)
	melt_table.rename(columns={
		melt_table.columns[0]:'sample',
		melt_table.columns[1]:'taxa',
		melt_table.columns[2]:'counts'}, inplace=True)

	# merge metadata and abundance data
	return pd.merge(metadata, melt_table, left_index=True, right_on='sample')


def _validate_columns(
	table: pd.DataFrame, 
	taxa_column: str,
	subject_column: str,
	pivot_column: str,
) -> None:
	if 'correlation_column' in table.columns:
		raise Exception("table cannot already has 'correlation_column'")
	if taxa_column not in table.columns:
		raise Exception('table does not have {c}'.format(c=taxa_column))
	if subject_column not in table.columns:
		raise Exception('table does not have {c}'.format(c=subject_column))
	if pivot_column not in table.columns:
		raise Exception('table does not have {c}'.format(c=pivot_column))


def correlation_plot(
	abundance_table: pd.DataFrame,
	metadata: pd.DataFrame,
	taxa_column: str,
	subject_column: str,
	pivot_column: str,
	base_value: str,
	save_loc: str,
) -> None:
	table = _melt_table(abundance_table, metadata)

	_validate_columns(
		table, taxa_column, subject_column, pivot_column)
	table['correlation_column'] = table[taxa_column] + table[subject_column]
	pivot_table = table.pivot(
		index='correlation_column', columns=pivot_column, values='counts')
	correlation_table = pivot_table.corr()
	correlation_table = correlation_table[[base_value]].reset_index()
	correlation_table = correlation_table[
		correlation_table[pivot_column] != base_value]

	# add error bars
	pivot_table['subject'] = pivot_table.index.str[-2:]
	pearsons = []
	for subject in pivot_table.subject.unique():
		data = pivot_table[pivot_table['subject'] == subject].corr()
		data['subject'] = subject
		pearsons.append(data)
	pearsons_table = pd.concat(pearsons)
	pearsons_table.reset_index(inplace=True)
	pearsons_table = pearsons_table[[pivot_column, base_value]]
	pearsons_table = pearsons_table[pearsons_table[pivot_column] != base_value]

	sns.set_context('talk')
	fig, ax = plt.subplots()
	fig.set_size_inches(5,5)
	sns.barplot(x=pivot_column, y=base_value, data=pearsons_table, ax=ax,
	            errwidth=2, capsize=.2, edgecolor='black')
	ax.set_ylabel('Pearson correlation to frozen swab', fontsize=16)
	ax.set_ylim(0,1)
	ax.set_xlabel('collection protocol')
	ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
	fig.savefig(save_loc)

	return correlation_table

meta = pd.read_csv('from-lisa-marotz/oral-16S-analysis/q2_saliva_metadata.tsv', sep='\t', index_col=0)
at = pd.read_csv('from-lisa-marotz/oral-16S-analysis/table-noMitoChlor-saliva-p6_.tsv', sep='\t', index_col=0)
correlation_plot(at, meta, 'taxa', 'host_subject_id', 'volume_ml', 'swbFroz', 'correlation_fig.png')
