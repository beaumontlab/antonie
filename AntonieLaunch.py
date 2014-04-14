import wx
import sys
import traceback
import  wx.lib.filebrowsebutton as filebrowse

class Frame(wx.Frame):
	def __init__(self, title):
		try: 
			wx.Frame.__init__(self, None, title=title, size=(850,650))
			
			self.fastq1 = filebrowse.FileBrowseButton(self, -1, size=(800, -1), changeCallback = self.fbbCallback, labelText= "FASTQ1 file", fileMask="FASTQ Files (*.fastq; *.fastq.gz)|*.fastq; *.gz|All files |*.*")
			self.fastq2 = filebrowse.FileBrowseButton(self, -1, size=(800, -1), changeCallback = self.fbbCallback, labelText= "FASTQ2 file", fileMask="FASTQ Files (*.fastq; *.fastq.gz)|*.fastq; *.gz|All files |*.*")
			self.fna = filebrowse.FileBrowseButton(	self, -1, size=(800, -1), changeCallback = self.fbbCallback, labelText="FASTA file", fileMask="FASTA Files (*.ffa; *.fna; *.fasta)|*.ffa;*.fna;*.fasta| All files |*.*")
			self.gff = filebrowse.FileBrowseButton(	self, -1, size=(800, -1), changeCallback = self.fbbCallback, labelText="GFF file", fileMask="GFF3 Files (*.gff; *.gff3)|*.gff;*.gff3| All files |*.*")

			
			self.process = None;
			self.Bind(wx.EVT_IDLE, self.OnIdle)			
			self.Bind(wx.EVT_END_PROCESS, self.OnProcessEnded)

			self.bam = wx.CheckBox(self, -1, "Create BAM file")
			self.unfound = wx.CheckBox(self, -1, "Create unfound.fastq")
			self.phix = wx.CheckBox(self, -1, "Exclude phiX")
			
			sizer = wx.BoxSizer(wx.VERTICAL)
			sizer.Add(self.fastq1, 0, wx.ALL, 5)
			sizer.Add(self.fastq2, 0, wx.ALL, 5)
			sizer.Add(self.fna, 0, wx.ALL, 5)
			sizer.Add(self.gff, 0, wx.ALL, 5)
			
			
			sizerh = wx.BoxSizer(wx.HORIZONTAL)
			sizerh.Add(self.bam, 0, wx.ALL, 5)
			sizerh.Add(self.unfound, 0, wx.ALL, 5)
			sizerh.Add(self.phix, 0, wx.ALL, 5)
			
			sizer.Add(sizerh, 0, wx.ALL, 5)
			
			self.goButton = wx.Button(self, -1, "Go")
			self.Bind(wx.EVT_BUTTON, self.OnGoBtn, self.goButton)
			self.goButton.Enable(False)

			sizer.Add(self.goButton, 0, wx.ALL, 5)
			self.out = wx.TextCtrl(self, -1, '', size=(800,300), style=wx.TE_MULTILINE|wx.TE_READONLY|wx.TE_RICH2)
			# Set a monospace font 
			tcFont = self.out.GetFont()
			fontPointSize = tcFont.GetPointSize() 
# No need to get the family ! 
			fontStyle     = tcFont.GetStyle() 
			fontWeight    = tcFont.GetWeight() 
			# To be thorough, [underline], [face] and [encoding] really should be  queried, too, 
			# but I'm lazy. 

			tcFont = wx.Font( fontPointSize, wx.FONTFAMILY_TELETYPE, fontStyle, 
					  fontWeight ) 
			self.out.SetFont( tcFont ) 

			sizer.Add(self.out, 0, wx.EXPAND|wx.ALL, 5);

			box = wx.BoxSizer()
			box.Add(sizer, 0, wx.ALL, 20)
			self.SetSizer(box)
		except:
			print "Unexpected error:", traceback.print_tb(sys.exc_info()[2])

	def fbbCallback(self, evt):
		if(self.fastq1.GetValue() != '' and self.fastq2.GetValue() != '' and
		   self.fna.GetValue() != '' and self.gff.GetValue() != ''):
			self.goButton.Enable(True)
		else:
			self.goButton.Enable(False)


	def OnGoBtn(self, evt):
		try:
			cmd = (".\\antonie.exe -1 \"" + self.fastq1.GetValue() + "\" -2 \"" +self.fastq2.GetValue() +"\" -r \""+self.fna.GetValue()+"\" -a \""+self.gff.GetValue()+"\"")
#			print "Go Button pressed, commandline: "+cmd
			
			self.process = wx.Process(self)
			self.process.Redirect();
			pid = wx.Execute(cmd, wx.EXEC_ASYNC, self.process)
		except:
			print "Unexpected error:", traceback.print_tb(sys.exc_info()[2])


	def OnIdle(self, evt):
		if self.process is not None:
			stream = self.process.GetInputStream()

			if stream.CanRead():
				text = stream.read()
				self.out.AppendText(text)

			stream = self.process.GetErrorStream()

			if stream.CanRead():
				text = stream.read()
				self.out.AppendText(text)


	def OnProcessEnded(self, evt):
		stream = self.process.GetInputStream()
		
		if stream.CanRead():
			text = stream.read()
			self.out.AppendText(text)

		stream = self.process.GetErrorStream()
		
		if stream.CanRead():
			text = stream.read()
			self.out.AppendText(text)

			
		self.process.Destroy()
		self.process = None
		self.goButton.Enable(True)



app = wx.App(redirect=True)
top = Frame("Antonie")
top.Show()
app.MainLoop()
