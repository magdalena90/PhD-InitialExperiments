
# 1. CREATE AN ENVIRONMENT WHERE TO INSTALL SCRAPY AND RUN CODE:

# First time:
# pip install --user pipenv

# PYTHON_BIN_PATH="$(python3 -m site --user-base)/bin"
# PATH="$PATH:$PYTHON_BIN_PATH"

# Activate environment and install scrapy:
# pipenv shell
# pipenv install Scrapy


# 2. CREATE PROJECT:
# scrapy startproject name-of-project (here it-s SFARI_scraping)
# Copy this file to name-of-project/name-of-project/spiders

# 3. MOVE INTO THE SCRAPY FOLDER
# cd name-of-project (here it-s SFARI_scraping)

# 4. RUN SPIDER:
# scrapy crawl name-assigned-to-genesInfoSpider (here it's genes_info)


# EXTRA: TO RUN SCRAPY ON SHELL:
# scrapy shell 'url_address'


import scrapy
import os
import csv
import re

class genesInfoSpider(scrapy.Spider):
	name = 'genes_info'
	start_urls = []

	with open('./../../../Data/SFARI/SFARI_genes_01-15-2019.csv') as csv_file:
		csv_reader = csv.reader(csv_file, delimiter=',')
		first_row = True
		for row in csv_reader:
			if(first_row):
				first_row = False
			else:
				start_urls.append('https://gene.sfari.org/database/human-gene/' + row[1])
	csv_file.close()

	def parse(self, response):

		# Create folder where to put output
		if not os.path.exists('genes'):
		    os.makedirs('genes')

		# Name output file
		gene = response.url.split("/")[-1]
		filename = 'genes/%s.csv' % gene

		# Write to file
		with open(filename, 'wb') as csv_file:

			writer = csv.writer(csv_file, delimiter='\t')

			whole_tbl = response.css('table')[0]

			# Header
			# header_tbl = whole_tbl.css('thead').css('th::text').getall()
			header_tbl = ['#','Type','Title','Author, year','pubmed ID', 'Autism Report', 'Associated Disorders']
			writer.writerow(header_tbl)

			# Body
			body_tbl = whole_tbl.css('tbody').css('tr')
			for row in body_tbl:
				row_content = row.css('td')
				n = row_content[0].css('td::text').extract()[0]
				Type = row_content[1].css('td::text').extract()[0]
				Title = row_content[2].css('td::text').extract()[0]
				Author = row_content[3].css('a::text').extract()[0]
				Author = re.sub('\t|\r|\n', '', Author)
				pubmedID = row_content[3].css('a::attr(href)').extract()[0].replace('https://www.ncbi.nlm.nih.gov/pubmed/','')
				Autism = row_content[4].css('td::text').extract()[0]
				Associated = row_content[5].css('td::text').extract()[0]

				new_row = [n, Type, Title, Author, pubmedID, Autism, Associated]

				writer.writerow(new_row)

		csv_file.close()

