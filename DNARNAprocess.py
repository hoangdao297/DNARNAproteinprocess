
class DNARNAprocess:
    aminoacidbase={'phe':['uuu','uuc'], 'leu':['uua','uug','cuu','cuc','cua','cug'],
        'ile':['auu','auc','aua'], 'val':['guu','guc','gua','gug'],
        'ala':['gcu','gcc','gca','gcg'],'thr':['acu','acc','aca','acg'],
        'pro':['ccu','ccc','cca','ccg'], 'ser':['ucu','ucc','uca','ucg','agu','agc'],
        'tyr':['uau','uac'], 'his':['cau','cac'], 'gln':['caa','cag'], 'asn':['aau','aac'],
        'lys':['aaa','aag'], 'asp':['gau','gac'], 'glu':['gaa','gag'], 'cys':['ugu','ugc'],
        'trp':['ugg'], 'arg':['cgu','cgc','cga','cgg','aga','agg'],
        'gly':['ggu','ggc','gga','ggg'], 'met':['aug']}
    shortcode={'A':'ala','C':'cys','D':'asp','E':'glu','F':'phe','G':'gly','H':'his','I':'ile',
        'K':'lys','L':'leu','M':'met','N':'asn','O':'pyl','P':'pro','Q':'gln','R':'arg','S':'ser',
        'T':'thr','U':'sec','V':'val','W':'trp','Y':'tyr'}
    stop=['uga','uaa','uag']
    start='aug'
    proteinlist=[]
    proteinlist2=[]
    codemain=""
    codeop=""
    def __init__(self,string):
        self.string=string
        self.type=self.dnaorrna() #0 is dna, 1 is rna
        self.opwaystr=self.processoppositeway()
        self.findingproteinformainstrand()
        self.findingproteinforopstrand()
        self.aminoacidmain=self.proteintoaminoacidmain()
        self.aminoacidop=self.proteintoaminoacidop()
        self.codemain=self.aminoacidcode(self.aminoacidmain)
        self.codeop=self.aminoacidcode(self.aminoacidop)
    
    def returnaminoacidmain(self):
        return self.aminoacidmain
    
    def returnaminoacidop(self):
        return self.aminoacidop

    def returncodemain(self):
        return self.codemain

    def returncodeop(self):
        return self.codeop
    
    def returnaminoacidcountmain(self):
        count=0
        for i in self.codemain:
            if (i!='-'): count+=1
        return count

    def returnaminoacidcountop(self):
        count=0
        for i in self.codeop:
            if (i!='-'): count+=1
        return count

    def dnaorrna(self):
        for i in self.string:
            if (i=='t'): return 0
            if (i=='u'): return 1
        return 0

    def torna(self,string):
        if (self.type==1): return string
        rna=''
        for i in string:
            if (i=='t'): rna+='u'
            else: rna+=i
        return rna

    def dnacomplement(self):
        if (self.type==1): return self.rnacomplement()
        compstr=''
        for i in self.string:
            if (i=='a'): compstr+='t'
            elif (i=='t'):compstr+='a'
            elif (i=='c'):compstr+='g'
            else:compstr+='c'
        return compstr

    def rnacomplement(self):
        if (self.type==0): return self.dnacomplement()
        compstr=''
        for i in self.string:
            if (i=='a'):compstr+='u'
            elif (i=='u'):compstr+='a'
            elif (i=='c'):compstr+='g'
            else:compstr+='c'
        return compstr

    def reversestr(self,string):
        reverse=''
        for i in range(len(string)):
            reverse+=string[len(string)-1-i]
        return reverse

    def processoppositeway(self):
        complementstr=self.rnacomplement()
        opstr=self.reversestr(complementstr)
        return opstr

    def findingproteinformainstrand(self):
        processstr=self.torna(self.string)
        #print(processstr)
        startcount=False
        i=2; proteintemp=""
        while (i<len(processstr)):
            tempstr=processstr[i-2]+processstr[i-1]+processstr[i]
            if (tempstr==self.start):
                if (len(self.proteinlist) ==0): self.proteinlist.append(self.start)
                if (len(self.proteinlist) >0 and self.proteinlist[len(self.proteinlist)-1]!=self.start): self.proteinlist.append(self.start)
                startcount=True;proteintemp="";i+=3;continue
            if (startcount==True):
                if (tempstr in self.stop):
                    self.proteinlist.append(proteintemp)
                    startcount=False
                else: 
                    proteintemp+=tempstr
                    i+=3; continue
            i+=1
        while (len(self.proteinlist)>0 and (self.proteinlist[len(self.proteinlist)-1]==self.start or self.proteinlist[len(self.proteinlist)-1]=='') ):
            self.proteinlist.pop(len(self.proteinlist)-1)
        return

    def findingproteinforopstrand(self):
        processstr=self.torna(self.opwaystr)
        startcount=False
        i=2; proteintemp=""
        while (i<len(processstr)):
            tempstr=processstr[i-2]+processstr[i-1]+processstr[i]
            if (tempstr==self.start):
                if (len(self.proteinlist2) ==0): self.proteinlist2.append(self.start)
                if (len(self.proteinlist2) >0 and self.proteinlist2[len(self.proteinlist2)-1]!=self.start): self.proteinlist2.append(self.start)
                startcount=True;proteintemp="";i+=3;continue
            if (startcount==True):
                if (tempstr in self.stop):
                    self.proteinlist2.append(proteintemp)
                    startcount=False
                else: 
                    proteintemp+=tempstr
                    i+=3; continue
            i+=1
        while (len(self.proteinlist2)>0 and (self.proteinlist2[len(self.proteinlist2)-1]==self.start or self.proteinlist2[len(self.proteinlist2)-1]=='')):
            self.proteinlist2.pop(len(self.proteinlist2)-1)
        return

    def proteintoaminoacidmain(self):
        protein= "".join(self.proteinlist)
        #print(protein)
        aminoacid=[]
        for i in range(0,len(protein),3):
            for aminoacidname, seq in self.aminoacidbase.items():
                if (protein[i:i+3] in seq):
                    aminoacid.append(aminoacidname); break
        #print(aminoacid)
        return aminoacid

    def proteintoaminoacidop(self):
        protein="".join(self.proteinlist2)
        aminoacid=[]
        for i in range(0,len(protein),3):
            for aminoacidname, seq in self.aminoacidbase.items():
                if (protein[i:i+3] in seq):
                    aminoacid.append(aminoacidname); break
        return aminoacid

    def aminoacidcode(self,aminoacidlist):
        code=""
        for i in aminoacidlist:
            for scode,name in self.shortcode.items():
                if (i==name):
                    code+=scode+'-';break
        #print(code)
        return code

def readandwritefile(filename):
    openfile='Input/'+filename+".txt"
    file=open(openfile,'r')
    lines=file.readlines()
    file.close()
    filewrite=open('Output/'+filename+'_result.txt','w')
    for i in range(len(lines)):
        string=validatestring(lines[i])
        result=DNARNAprocess(string)
        filewrite.write("Aminoacid in main strand: ")
        filewrite.write("-".join(result.returnaminoacidmain()))
        filewrite.write("\n")
        filewrite.write("Aminoacid code in main strand: "+result.returncodemain()+"\n")
        filewrite.write("Amino acid count in main strand: "+str(result.returnaminoacidcountmain())+'\n')
        filewrite.write("Aminoacid in op strand: ")
        filewrite.write("-".join(result.returnaminoacidop()))
        filewrite.write("\n")
        filewrite.write("Aminoacid code in op strand: "+result.returncodeop()+"\n")
        filewrite.write("Amino acid count in op strand: "+str(result.returnaminoacidcountop())+'\n\n')
    filewrite.close()

def validatestring(string):
    return string.lower()
    
readandwritefile('input2')
#rna='augaugcccuagccaugcccuagaugaugauguag'
#process=DNARNAprocess(rna)
#print(process.returncodemain())
#print(process.returnaminoacidcountmain())

