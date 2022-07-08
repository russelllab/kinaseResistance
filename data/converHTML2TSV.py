import os

num = 0; l = ''; row = []
for line in open('kinases.html', 'r'):
	if '<td>' in line and '</td>' in line:
		row.append(str(line.replace(' ','').replace('\n', '').replace('<td>', '').replace('</td>', '')))
		num += 1
		if num == 8:
			print ('\t'.join(row))
			l += '\t'.join(row) + '\n'
			row = []
			num = 0
			
print (l)
open('kinases.tsv', 'w').write(l)
