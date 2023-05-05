import pandas as pd
	
def query_nea(col='*', fmt='csv', tab='ps', const=None):
	"""
		Query the NEA database with the TAP protocol.
		
		Inputs
		------
			
			col: str
				Columns to select from the NEA database. Default='*' for all.
				
			fmt: str
				Output format of the results. Default='csv'.
				
			tab: str
				NEA table to query. Default='ps'
				
			const: dict
				Constraints to be applied when querying. This dict should
				contain 3 key-value pairs:
				
					'col':[]
					'oper':[]
					'val':[]
				
				where 'col' is the column to constrain, 'oper' is the relational
				operator, and 'val' is the criterion. Default=None.
				
		Outputs
		-------
		
			df: pandas.DataFrame
				Object containing the result of the NEA query.
				
	"""
	
	# Construct the url
	url = construct_nea_url(col, fmt, tab, const=const)
	
	try:
		df = pd.read_csv(url)
	except Exception as ee:
		print(ee)
		print('The URL sent was: {}'.format(url))
		
	return df	

def construct_nea_url(col, fmt, tab, const=None):
	"""
		Construct the URL for the call to the NEA database.
		
		Inputs
		------
		
			col: str
				Columns to select from the NEA database. 
				
			fmt: str
				Output format of the results. 
				
			tab: str
				NEA table to query. 
				
			const: dict
				Constraints to be applied when querying. This dict should
				contain 3 key-value pairs:
				
					'col':[]
					'oper':[]
					'val':[]
				
				where 'col' is the column to constrain, 'oper' is the relational
				operator, and 'val' is the criterion. Default=None.
				
		Outputs
		-------
		
			url: str
				Concatenated URL.
	"""
	
	base_url = 'https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query=select+'
	url = base_url + col
	url += '+from+'
	url += tab
	if const is not None:
		url += '+where+'
		for i in range(len(const['col'])):
			url += const['col'][i] + '+'
			url += const['oper'][i] + '+'
			url += const['val'][i]
			if i < (len(const['col']) - 1):
				url += '+and+'
	url += '&format={}'.format(fmt)

	return url