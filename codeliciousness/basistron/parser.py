from html.parser import HTMLParser

class TableParser(HTMLParser):

	def __init__(self):
		super().__init__()
		self.td = False
		self.th = False
		self.skip = False
		self.live = False
		self.current_table = []
		self.current_row = []
		self.current_cell = []
		self.tables = []

	def handle_starttag(self, tag, attrs):
		if tag == "td":
			self.td = True
			self.live = True
		if tag == "th":
			self.th = True
			if len(attrs) and attrs[0][0] == "rowspan":
				self.skip = True
		if tag == "a" and self.td:
			self.current_cell.append(attrs[0][1])
    
	def handle_data(self, data):
		if not self.skip and (self.td or self.th):
			self.current_cell.append(data.strip())

	def handle_endtag(self, tag):
		if tag == "td":
			self.td = False
		if tag == "th":
			self.th = False
		if tag in ["td", "th"]:
			if not self.skip:
				self.current_row.append(self.current_cell[::-1])
			self.current_cell = []
			self.skip = False
		if tag == "tr":
			self.current_table.append(self.current_row)
			self.current_row = []
		if tag == "table":
			self.tables.append(self.current_table)
			self.current_table = []


