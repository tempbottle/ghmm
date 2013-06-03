#!/usr/bin/env python
################################################################################
#
#       This file is part of the General Hidden Markov Model Library,
#       GHMM version __VERSION__, see http://ghmm.org
#
#       file:    xmlutil.py
#       authors: Wasinee Rungsarityotin, Janne Grunau
#
#       Copyright (C) 1998-2004 Alexander Schliep
#       Copyright (C) 1998-2001 ZAIK/ZPR, Universitaet zu Koeln
#       Copyright (C) 2002-2004 Max-Planck-Institut fuer Molekulare Genetik,
#                               Berlin
#
#       Contact: schliep@ghmm.org
#
#       This library is free software; you can redistribute it and/or
#       modify it under the terms of the GNU Library General Public
#       License as published by the Free Software Foundation; either
#       version 2 of the License, or (at your option) any later version.
#
#       This library is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#       Library General Public License for more details.
#
#       You should have received a copy of the GNU Library General Public
#       License along with this library; if not, write to the Free
#       Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
#
#
################################################################################


from DataStructures import Point2D,EdgeWeight
from Graph import Graph, SubGraph

from xml.dom.minidom import *
import EditObjectAttributesDialog
from EditObjectAttributesDialog import ValidatingString, ValidatingInt, ValidatingFloat, PopupableInt, Probability, DefaultedInt, DefaultedString

import string
import types
import copy

import sys 

import logging
log = logging.getLogger("xmlutil.py")

def typed_assign(var, val):
    result = type(var)(val)
    result.__dict__ = var.__dict__
    #result.__dict__ = copy.copy(var.__dict__)
    return result

def listFromCSV(s, type):
    return map(type,string.split(s,','))

def csvFromList(list, perRow = None):
    if perRow == None:
        return string.join(map(str,list), ', ')
    else:
        result = ""
        for start in xrange(0, len(list), perRow):
            result += string.join(map(str,list[start:start+perRow]), ', ') + ',\n'
        return result[0:len(result)-2]

def writeContents(XMLDoc, XMLNode, data):
    contents = XMLDoc.createTextNode("%s" % data)
    XMLNode.appendChild(contents)

def writeData(XMLDoc, XMLNode, dataKey, dataValue):
    data = XMLDoc.createElement("data")
    data.setAttribute('key', "%s" % dataKey)
    contents = XMLDoc.createTextNode("%s" % dataValue)
    data.appendChild(contents)
    XMLNode.appendChild(data)

def writeXMLData(XMLDoc, XMLNode, dataKey, XMLData):
    data = XMLDoc.createElement("data")
    data.setAttribute('key', "%s" % dataKey)
    data.appendChild(XMLData)
    XMLNode.appendChild(data)

def writeXMLTextNode(XMLDoc, XMLNode, keyName, XMLData):
    data = XMLDoc.createElement(keyName)
    contents = XMLDoc.createTextNode("%s" % XMLData)
    data.appendChild(contents)
    XMLNode.appendChild(data)

class NamedDistributions:
    def __init__(self, itsHMM):
        self.initialize()
        self.itsHMM = itsHMM
        self.code2name = {-1:'None'}
        self.name2code = {'None':-1}
        self.maxCode = 0

    def initialize(self):
        self.dist = {}
        self.order = {}
       
    def addDistribution(self, name, order, p): 
        self.dist[name] = p
        self.order[name] = order
        self.code2name[self.maxCode] = name
        self.name2code[name] = self.maxCode
        self.maxCode += 1
        # print self.dist
        
    def deleteDistribution(self, name):
        del self.dist[name]
        del self.order[name]
        del self.code2name[self.name2code[name]]
        del self.name2code[name]
   
    def fromDOM(self, XMLNode):
        self.initialize()
        datas = XMLNode.getElementsByTagName("hmm:background")
        for data in datas:
            dataKey = data.attributes['key'].nodeValue
            dataOrder = int(data.attributes['order'].nodeValue)
            dataValue = ""
            for child in data.childNodes:
                dataValue += child.nodeValue
            p = listFromCSV(dataValue, types.FloatType)
            self.addDistribution(dataKey, dataOrder, p)

    def toDOM(self, XMLDoc, XMLNode):
        for name in self.dist.keys():
            #print "background: name = ", name, self.dist[name], self.order[name]
            background_elem = XMLDoc.createElement("hmm:background")
            background_elem.setAttribute('key', "%s" % name)
            background_elem.setAttribute('order', "%s" % self.order[name])
            if self.order[name] == 0:
                contents = XMLDoc.createTextNode(csvFromList(self.dist[name]))
            else:
                contents = XMLDoc.createTextNode(csvFromList(self.dist[name],
                                                             self.itsHMM.hmmAlphabets[0].size()))               
            background_elem.appendChild(contents)
            XMLNode.appendChild(background_elem)

    def names(self):
        return self.dist.keys()

   

class XMLElementWriter:
    
    def __init__(self):
        import string
        self._string = string
        self.element_keys = ['id', 'hmm:low', 'hmm:high', 'hmm:type', 'hmm:id','hmm:alphabet_id', 'hmm:name', 'xmlns', 'xmlns:gd', 'xmlns:hmm', 'for', 'gd:type', 'code', 'key', 'order', 'x', 'y', 'source', 'target', 'name', 'domain']
        
    def _write_data(self, writer, data):
        "Writes datachars to writer."
        replace = self._string.replace
        data = str(data)
        data = replace(data, "&", "&amp;")
        data = replace(data, "<", "&lt;")
        data = replace(data, "\"", "&quot;")
        data = replace(data, ">", "&gt;")
        writer.write(data)

    def writexml(self, XMLNode, writer, indent="", addindent="", newl=""): 
        """This version of write xml makes sure text nodes are 
        surrounded by tags for easier reading rather than being 
        on lines by themselves."""  

        # indent = current indentation
        # addindent = indentation to add to higher levels
        # newl = newline string
        
        if ( XMLNode.nodeType == XMLNode.TEXT_NODE ): 
            self._write_data(writer, "%s%s%s"%(indent, XMLNode.data, newl))
        # also handle comments
        elif (XMLNode.nodeType == XMLNode.COMMENT_NODE):
            writer.write("%s<!--" % indent)
            self._write_data(writer, "%s" % XMLNode.nodeValue)
            writer.write("-->%s" % newl)
        else:
            writer.write(indent+"<" + XMLNode.nodeName)
            
            # build attribute list
            a_names = []
            try:
                for key in self.element_keys:
                    if ( XMLNode.getAttribute(key) != ""):
                        a_names.append( key )
                        a_names.sort()
            except AttributeError:
                a_names = []
        
            for a_name in a_names:
                writer.write(" %s=\"" % a_name)
                self._write_data(writer, XMLNode.getAttribute(a_name))
                writer.write("\"")
            if XMLNode.childNodes:
                writer.write(">%s"%(newl))
                for node in XMLNode.childNodes:
                    if node.nodeType!=node.TEXT_NODE: 
                        self.writexml(node,writer,indent+addindent,addindent,newl)
                    else:
                        writer.seek(writer.tell()-1)
                        self.writexml(node,writer,"",addindent,"")

                if XMLNode.childNodes[-1].nodeType!= XMLNode.TEXT_NODE:
                    writer.write("%s</%s>%s" % (indent,XMLNode.nodeName,newl))
                else:
                    writer.write("</%s>%s" % (XMLNode.nodeName,newl))
            else:
                writer.write("/>%s"%(newl))
                     

#
# Notice for this function:
# minidom.Element.toprettyxml is not so pretty and its output format is not compatible
# with XMLIO parser. The particular problem with minidom.toprettyxml is that
# it put the text data of a text node on a new line, instead of immediately after the element tag.
# Because XMLIO cannot parse this format, thus we need our own pretty print program 
#
def toprettyxml( XMLDoc ):
    # we can't use cStringIO since it doesn't support Unicode strings
    from StringIO import StringIO
    writer = StringIO()
    prettydoc = XMLElementWriter()
    writer.write('<?xml version="1.0" ?>\n')
    # remove empty text nodes to make XMLElementWriter work correctly
    removeEmptyTextNodes(XMLDoc)
    for node in XMLDoc.childNodes:
        prettydoc.writexml(node, writer, "","  ", "\n")
    return writer.getvalue();

def removeEmptyTextNodes(XMLNode):
    # some minidom implementations include lots of whitespace text nodes
    # they prevent the XMLWriter from writing nice XML
    if (XMLNode.nodeType == XMLNode.TEXT_NODE):
        if (XMLNode.nodeValue.strip() == ""):
            XMLNode.parentNode.removeChild(XMLNode)
    else:
        # have to copy the list because the removeChild() modifies the list
        import copy
        childNodes = copy.copy(XMLNode.childNodes)
        for child in childNodes:
            removeEmptyTextNodes(child)

class DOM_Map:
    def __init__(self):
        self.initialize()

    def initialize(self):
        self.name = {}
        self.desc = {}
        self.hasDesc = None
        self.name2code = {}
        
    def addCode(self, code, name, desc = None):
        self.name[code] = name
        if desc != None:
            self.desc[code] = desc
            self.hasDesc = 1
        self.name2code[name] = code

    def low(self):
        if len(self.name.keys()) > 0:
            return min(self.name.keys())
        else:
            return 0
                
    def high(self):
        if len(self.name.keys()) > 0:
            return max(self.name.keys())
        else:
            return 0
    
    def fromDOM(self, XMLNode):
        pass

    def symbolsFromDom(self, XMLNode):
        symbols = XMLNode.getElementsByTagName("symbol")
        
        for symbol in symbols:
            symbolCode = ValidatingInt(int(symbol.getAttribute("code")))
            symbolName = ValidatingString(symbol.firstChild.nodeValue)
            symbolDesc = symbol.getAttribute("desc")
            if symbolDesc != None:
                self.addCode(symbolCode, symbolName, ValidatingString(symbolDesc))
            else:
                self.addCode(symbolCode, symbolName)
                
    def toDOM(self, XMLDoc, XMLNode):
        XMLNode.setAttribute('hmm:low', "%s" % self.low())
        XMLNode.setAttribute('hmm:high', "%s" % self.high())
        map = XMLDoc.createElement("map")  
        for key in self.name.keys():
            symbol = XMLDoc.createElement("symbol")
            symbol.setAttribute('code', "%s" % key)
            if self.hasDesc and self.desc[key] != "":
                symbol.setAttribute('desc', "%s" % self.desc[key])
            writeContents(XMLDoc, symbol, "%s" % self.name[key])
            map.appendChild(symbol)
        XMLNode.appendChild(map)
   
    def buildList(self):
	return self.name.keys()
     
# -------------------------------------------
#  Exceptions

class HMMEdError(Exception):
    def __init__(self, message):
	print "\n\n Unknown error types. Please report \n\n"
	
class NotValidHMMType(HMMEdError):
    def __init__(self,message):
       print "\n\n Probabilities missing xception: " + str(message) + "\n"

class AlphabetErrorType(HMMEdError):
    def __init__(self,message):   
        print "\n\n Alphabet exception: " + str(message) + "\n"
	
	
class DiscreteHMMAlphabet(DOM_Map):
    def __init__(self):
        DOM_Map.__init__(self)
        self.hmm_type = 'discrete'
        self.id = 0

    def fromDOM(self, XMLNode):
        """Take dom subtree representing a <hmm:alphabet</hmm:alphabet> element"""
        self.initialize()
        # Not reading: hmm:low hmm:high
        if XMLNode.getAttribute("hmm:type") == self.hmm_type:
            self.symbolsFromDom(XMLNode)
            id = XMLNode.getAttribute("hmm:id")
            if (id != ""):
                self.id = int(id)
        else:
            print "DiscreteHMMAlphabet wrong type %s" % XMLNode.getAttribute("hmm:type") 

    def toDOM(self, XMLDoc, XMLNode):
        hmmalphabet = XMLDoc.createElement("hmm:alphabet")
        hmmalphabet.setAttribute('hmm:type', 'discrete')
	hmmalphabet.setAttribute('hmm:low', "%s" % self.low())
        hmmalphabet.setAttribute('hmm:high', "%s" % self.high())
        hmmalphabet.setAttribute('hmm:id', str(self.id))
        map = XMLDoc.createElement("map")  
        for key in self.name.keys():
            symbol = XMLDoc.createElement("symbol")
            symbol.setAttribute('code', "%s" % key)
            if self.hasDesc and self.desc[key] != "":
                symbol.setAttribute('desc', "%s" % self.desc[key])
            writeContents(XMLDoc, symbol, "%s" % self.name[key])
            map.appendChild(symbol)
        hmmalphabet.appendChild(map)
	#  DOM_Map.toDOM(self, XMLDoc, hmmalphabet)
        XMLNode.appendChild(hmmalphabet)

    def toGHMM(self, XMLDoc, XMLNode):
        hmmalphabet = XMLDoc.createElement("alphabet")
        for key in self.name.keys():
            alphabet = XMLDoc.createElement("symbol")
            alphabet.setAttribute('id', "%s" % key)
            hmmalphabet.appendChild(alphabet)
        XMLNode.appendChild(hmmalphabet)
      
   
    def size(self):
        return len(self.name.keys())
    
    def buildList(self):
	return self.name.keys()

    def buildAlphabets(self, nrOfSymbols):
	""" Only fixed to 5 symbols for the moment """ 
	alphas = range(1, nrOfSymbols+1)
	alphas = map( lambda x: 'a'+ str(x), alphas)
	for code in range(1,nrOfSymbols+1):
	    self.addCode( code, alphas[code-1], desc = None)
 


class HMMClass(DOM_Map):
    def __init__(self):
        DOM_Map.__init__(self)
        self.code2name = {-1:'None'}
        self.name2code = {'None':-1}
        self.maxCode = 0
    
    def fromDOM(self, XMLNode):
        """Take dom subtree representing a <hmm:class></hmm:class> element"""
        self.initialize()
        self.symbolsFromDom(XMLNode)
        # copy self.name to self.code2name
        for key in self.name.keys():            
            self.code2name[key] = self.name[key]
        
    def toDOM(self, XMLDoc, XMLNode):        
        hmmclass = XMLDoc.createElement("hmm:class")   
        DOM_Map.toDOM(self, XMLDoc, hmmclass)
        XMLNode.appendChild(hmmclass)

    def size(self):
        return len(self.name.keys())

            
class HMMState:

    def __init__(self, nodeIndex, itsHMM):

        self.initial  = Probability("0.0")
        self.label    = ValidatingString("None")
        self.itsHMM   = itsHMM
	# print type(self.label)
	
        self.index = nodeIndex # The node index in the underlying graph
        self.id    = DefaultedInt() # identification by the canvas, not always the same
        self.state_class = PopupableInt(-1)
        self.state_class.setPopup(itsHMM.hmmClass.code2name, itsHMM.hmmClass.name2code, 10)
        # XXX self.state_class.setPopup(itsHMM.hmmClass.name, itsHMM.hmmClass.name2code, 10)

	self.order = DefaultedInt()
        self.order.setDefault(1, 0) # the default state order is 0

        self.emissions = []

        self.tiedto = DefaultedString()
        self.tiedto.setDefault(1, '')
        self.desc = self.id

        self.reading_frame = PopupableInt(-1)
        code2name = {-1:'None', 0:'0', 1:'1', 2:'2'}
        name2code = {'None':-1, '0':0, '1':1, '2':2}
        self.reading_frame.setPopup(code2name, name2code, 4)

        self.duration = DefaultedInt()
        self.duration.setDefault(1, 0)

        self.background = PopupableInt(-1)
        self.background.setPopup(self.itsHMM.backgroundDistributions.code2name, self.itsHMM.backgroundDistributions.name2code, 10)

        self.editableAttr = ['label', 'state_class', 'initial', 'order', 'background', 'offsetX', 'offsetY', 'alphabet_id']
        self.xmlAttr = self.editableAttr + ['ngeom', 'emissions']
        # pair HMM stuff
        self.alphabet_id = ValidatingInt(0)
        self.offsetX = ValidatingInt(1)
        self.offsetY = ValidatingInt(0)
        self.kclasses = ValidatingInt(1)
        self.transitionFunction = ValidatingInt(-1)
        
    editableAttr = ['label', 'state_class', 'initial', 'order', 'background', 'offsetX', 'offsetY', 'alphabet_id']
    xmlAttr = editableAttr + ['ngeom', 'emissions']
    
    # ['id', 'state_class', 'label', 'order', 'initial', 'tiedto', 'reading_frame', 'duration', 'background']

    #
    # Set_<Attribute> Methods: for integration with class editobj.Editor.set_value()
    # the name of the method must be in the form of set_<attr name>(self,value).
    # Otherwise, EditObj cannot propogate new values. (See the base method editobj.EntryEditor.set_value())
    # When necessary, we also update the Graph here.
    #
    def set_label(self, value):
        # Get the id2index from the Graph
        oldv = self.itsHMM.id2index[self.id]
        self.label = typed_assign(self.label, value)
        
        # We only show the label out of the editable items
        self.itsHMM.G.labeling[oldv] = "%s\n%s" % (self.id, self.label) # XXX Hack Aaaargh!
        self.itsHMM.itsEditor.UpdateVertexLabel(oldv, 0)


    def set_initial(self, value):
        self.initial = typed_assign(self.initial, value)

    def set_offsetX(self, value):
        self.offsetX = typed_assign(self.offsetX, value)

    def set_offsetY(self, value):
        self.offsetY = typed_assign(self.offsetY, value)

    def set_alphabet_id(self, value):
        self.alphabet_id = typed_assign(self.alphabet_id, value)

    def fromDOM(self, XMLNode):
        
        self.id = typed_assign(self.id, int(XMLNode.attributes['id'].nodeValue)) # state's id
        
        self.index = self.itsHMM.G.AddVertex()
        
        datas = XMLNode.getElementsByTagName("data")
        for data in datas:
            dataKey = data.attributes['key'].nodeValue
            dataValue = data.firstChild.nodeValue

            #print dataValue
            
            #  if dataValue == None: # use default Value
            #     self.state_class = typed_assign(self.state_class, int(0))
            #  else:
            #      self.state_class = typed_assign(self.state_class, int(dataValue))

            if dataKey == 'class':
                self.state_class = typed_assign(self.state_class, int(dataValue))
                    
            elif  dataKey == 'label':
                self.label = type(self.label)(dataValue.encode('ascii', 'replace'))

            elif  dataKey == 'order':
                if dataValue == None: # use default value
                    self.order = typed_assign(self.order, self.order.defaultValue)
                    self.order.useDefault = 1
                else:
                    self.order = typed_assign(self.order, int(dataValue))
                    self.order.useDefault = 0

            elif  dataKey == 'initial':
                self.initial = typed_assign(self.initial, float(dataValue))
                
            elif  dataKey == 'tiedto':
                
                if dataValue == None: # use default value
                    self.tiedto = typed_assign(self.tiedto, self.tiedto.defaultValue)
                    self.tiedto.useDefault = 1
                else:
                    self.tiedto = typed_assign(self.tiedto, dataValue.encode('ascii', 'replace'))
                    self.tiedto.useDefault = 0

            elif dataKey == 'reading-frame':
                self.reading_frame = typed_assign(self.reading_frame, int(dataValue))

            elif dataKey == 'background':
                self.background = typed_assign(self.background, self.itsHMM.backgroundDistributions.name2code[dataValue])

            elif dataKey == 'duration':
                self.duration = typed_assign(self.duration, int(dataValue))
                self.duration.useDefault = 0
                                    
            elif dataKey == 'ngeom':
                # We only use pos
                pos = XMLNode.getElementsByTagName('pos')[0] # Just one pos ...                
                self.pos = Point2D(float(pos.attributes['x'].nodeValue),
                                   float(pos.attributes['y'].nodeValue))
                
            elif dataKey == 'emissions':
                # collect all strings from childnodes
                dataValue = ""
                for child in data.childNodes:
                    dataValue += child.nodeValue
                self.emissions = listFromCSV(dataValue, types.FloatType)
                #print self.emissions

            elif dataKey == 'alphabet_id':
                self.alphabet_id = ValidatingInt(dataValue)

            elif dataKey == 'offset_x':
                self.offsetX = ValidatingInt(dataValue)

            elif dataKey == 'offset_y':
                self.offsetY = ValidatingInt(dataValue)

            elif dataKey == 'kclasses':
                self.kclasses = ValidatingInt(dataValue)

            elif dataKey == 'transitionfunction':
                self.transitionFunction = ValidatingInt(dataValue)
                    
            else:
                print "HMMState.fromDOM: unknown key %s of value %s" % (dataKey, dataValue)
        

    def toDOM(self, XMLDoc, XMLNode, initial_sum):
        node = XMLDoc.createElement("node")
        node.setAttribute('id', "%s" % self.id)

        # Mandatory elems
        writeData(XMLDoc, node, 'label', self.label)
        writeData(XMLDoc, node, 'class', self.state_class)
        writeData(XMLDoc, node, 'initial', self.initial / initial_sum)

        pos_elem = XMLDoc.createElement("pos")
        try:
            pos = self.itsHMM.G.embedding[self.index] 
            pos_elem.setAttribute('x', "%s" % pos.x)
            pos_elem.setAttribute('y', "%s" % pos.y)
            writeXMLData(XMLDoc, node, 'ngeom', pos_elem)
        except:
            pos_elem.setAttribute('x', "%s" % 100)
            pos_elem.setAttribute('y', "%s" % 50)
            writeXMLData(XMLDoc, node, 'ngeom', pos_elem)

        writeData(XMLDoc, node, 'offset_x', self.offsetX)

        writeData(XMLDoc, node, 'offset_y', self.offsetY)

        writeData(XMLDoc, node, 'alphabet_id', self.alphabet_id)

        writeData(XMLDoc, node, 'kclasses', self.kclasses)

        writeData(XMLDoc, node, 'transitionfunction', self.transitionFunction)
            
        if not self.order.useDefault:
            writeData(XMLDoc, node, 'order', self.order)

        if self.reading_frame != -1:
            writeData(XMLDoc, node, 'reading-frame', self.reading_frame)

        if self.background != -1:
            writeData(XMLDoc, node, 'background', self.itsHMM.backgroundDistributions.code2name[self.background])

        if not self.duration.useDefault:
            writeData(XMLDoc, node, 'duration', self.duration)
           
        if not self.tiedto == '':
            writeData(XMLDoc, node, 'tiedto', self.tiedto)
            self.emissions = self.itsHMM.state[self.itsHMM.id2index[int(self.tiedto)]].emissions # XXX            

        if self.order.useDefault:
            order = 0
        else:
            order = self.order

        # XXX Produce uniform emission probs, if we dont have the correct number of
        # parameters
            
        size = self.itsHMM.hmmAlphabets[self.alphabet_id].size()**(order+1)
        if (self.itsHMM.modelType == "pairHMM"):
            if (self.offsetX != 0 and self.offsetY != 0):
                size = size**2
        if len(self.emissions) != size:
            print "in state %s len(emissions) = %i size should be %i" % (self.id, len(self.emissions), size)
            tmp = [1.0/self.itsHMM.hmmAlphabets[self.alphabet_id].size()] * self.itsHMM.hmmAlphabet.size()
            if order == 0:
                self.emissions = tmp
            else:
                self.emissions = tmp * self.itsHMM.hmmAlphabets[self.alphabet_id].size()**order
                    
                
        if order > 0:
            writeData(XMLDoc, node, 'emissions', csvFromList(self.emissions,
                                                             self.itsHMM.hmmAlphabets[self.alphabet_id].size()))
        else:
            writeData(XMLDoc, node, 'emissions', csvFromList(self.emissions))
            
        XMLNode.appendChild(node)

    def toGHMM(self, XMLDoc, XMLNode, initial_sum):
        node = XMLDoc.createElement("state")
        node.setAttribute('id', "%s" % self.id)

        writeXMLTextNode(XMLDoc, node, 'initial', self.initial / initial_sum)
        # ignore order
        writeXMLTextNode(XMLDoc, node, 'emission', string.join(map(str,self.emissions),'\n'))
        XMLNode.appendChild(node)
        

    def fromDiscreteState( self, id, pi, B, label, order, tiedto=None, background=False ):
        """ Convert from ghmm.DiscreteEmissionHMM to the member attributes
        """
        self.id = id
        self.index = self.itsHMM.G.AddVertex()
        self.initial = typed_assign(self.initial, float(pi))
        self.label   = type(self.label)(label)	
	self.order = DefaultedInt(order)

        self.order = typed_assign(self.order, int(order))
        self.order.useDefault = 0
        self.emissions = B

        if tiedto is None:
            self.tiedto.useDefault = 1
        else:
            self.tiedto = typed_assign(self.tiedto, tiedto)
            self.tiedto.useDefault = 0

        if background is True:
            self.background = PopupableInt(-1)
            self.background.setPopup(self.itsHMM.backgroundDistributions.code2name, self.itsHMM.backgroundDistributions.name2code, 10)


class HMM:    
    def __init__(self, XMLFileName = None, G = None):
        # self.itsEditor = itsEditor
        if ( G is None ):
            self.G = Graph()
        else:
            self.G = G

        self.G.directed = 1
        self.G.euclidian = 0
        self.G.simple = 0
        self.Pi = {}
        self.id2index = {}
            
        # self.hmmAlphabet = DiscreteHMMAlphabet()
        self.hmmClass    = HMMClass()

        # in the case of pair HMMs we have several
        self.hmmAlphabets = {}
        self.transitionFunctions = {}
        
        self.editableAttr = {}
        self.editableAttr['HMM'] = ['desc']
        self.desc = ValidatingString()       

        self.state = {}

        self.modelType = 0
        self.name = "NoName"
        
        self.backgroundDistributions = NamedDistributions(self)

        self.DocumentName = "graphml"
        if XMLFileName != None:
            self.OpenXML(XMLFileName)

    def Clear(self):
        self.G.Clear()
        self.Pi = {}
        self.id2index = {}
            
        # self.hmmAlphabet = DiscreteHMMAlphabet()
        self.hmmAlphabets = {}
        self.hmmClass    = HMMClass()
        self.backgroundDistributions = NamedDistributions(self)
        
        self.editableAttr = {}
        self.editableAttr['HMM'] = ['desc']
        self.desc = ValidatingString()       
        self.state = {}
        self.DocumentName = "graphml"  

    def AddState(self, index, label='None'):
        state = HMMState(-1, self)
        if self.id2index.keys() != []:
            state.id = max(self.id2index.keys()) + 1
        else:
            state.id = 1
        state.index = index
        self.id2index[state.id] = state.index
        self.state[state.index] = state # XXX Use canvas id
        state.label = typed_assign(state.label, state.id)
        self.G.labeling[state.index] = "%s" % (state.label)
        return state.index
        
    def DeleteState(self, index):
	""" The method only deletes a map between index and its state object.
	    The caller must delete the corresponding vertex in the owner Graph self.G. """
	del self.id2index[self.state[index].id]
	del self.state[index]

    def fromDOM(self, XMLNode):
        
        # self.hmmClass.fromDOM(XMLNode.getElementsByTagName("hmm:class")[0]) 
        for tag in XMLNode.getElementsByTagName("hmm:class"):
            self.hmmClass.fromDOM(tag)

        nameNodes = XMLNode.getElementsByTagName("hmm:name")
        if (len(nameNodes) > 0):
            self.modelType = nameNodes[0].firstChild.nodeValue

        # model type node
        modelTypeNodes = XMLNode.getElementsByTagName("hmm:modeltype")
        if (len(modelTypeNodes) > 0):
            self.modelType = modelTypeNodes[0].firstChild.nodeValue
        if (self.modelType == "pairHMM"):
            alphabetNodes = XMLNode.getElementsByTagName("hmm:alphabet")
            for alphabetNode in alphabetNodes:
                alphabet = DiscreteHMMAlphabet()
                alphabet.fromDOM(alphabetNode)
                self.hmmAlphabets[alphabet.id] = alphabet
            transitionFunctionNodes = XMLNode.getElementsByTagName("hmm:transitionfunction")
            for transitionFunctionNode in transitionFunctionNodes:
                transitionFunction = TransitionFunction()
                transitionFunction.fromDom(transitionFunctionNode)
                self.transitionFunctions[transitionFunction.id] = transitionFunction
        else:
            # If it is no pair hmm One "hmm:alphabet" XML element
            self.hmmAlphabets[0] = DiscreteHMMAlphabet()
            self.hmmAlphabets[0].fromDOM(XMLNode.getElementsByTagName("hmm:alphabet")[0])

        self.backgroundDistributions.fromDOM(XMLNode)

        nodes = XMLNode.getElementsByTagName("node")
        for n in nodes:
            state = HMMState(-1, self)
            state.fromDOM(n)
            self.state[state.index] = state # key must be string
            self.id2index[state.id] = state.index
            self.G.embedding[state.index] = state.pos
            self.G.labeling[state.index] = "%s\n%s" % (state.id, state.label) # XXX Hack Aaaargh!

        edges = XMLNode.getElementsByTagName("edge")
        # nr_classes = int(self.hmmClass.high()-self.hmmClass.low())+1
        nr_classes = 1
        # search in all states for the maximal kclasses
        for s in self.state.values():
            if (s.kclasses > nr_classes):
                nr_classes = s.kclasses
                
        for i in range(nr_classes):
            self.G.edgeWeights[i] = EdgeWeight(self.G)

        for edge in edges:
            i = self.id2index[int(edge.attributes['source'].nodeValue)]
            j = self.id2index[int(edge.attributes['target'].nodeValue)]
            source = self.state[i]
            datas = edge.getElementsByTagName("data")
            for data in datas:
                dataKey = data.attributes['key'].nodeValue
                # dataValue = data.firstChild.nodeValue

            if dataKey == 'prob':
                #p = float(dataValue)
                # collect all strings from childnodes
                dataValue = ""
                for child in data.childNodes:
                    dataValue += child.nodeValue
                p = listFromCSV(dataValue, types.FloatType)
                self.G.AddEdge(i, j)
                if len(p) == 1: # only one class
                    for cl in range(source.kclasses - 1):
                        p.append(0.0)
                        
                for cl in range(source.kclasses):
                    self.G.edgeWeights[cl][(i,j)] = p[cl]

    def modelCheck(self):
	
        # Compute sums of initial probabilities for renormalization 
        initial_sum = 0.0
        for s in self.state:
            initial_sum = initial_sum + self.state[s].initial

	if initial_sum == 0.0:
	    raise NotValidHMMType("Initial state is not specified.")
	    
	if (len(self.hmmAlphabets) == 0):
	    raise AlphabetErrorType("Alphabet object is empty. You must create alphabet before saving.")
	
    def toDOM(self, XMLDoc, XMLNode):
        graphml = XMLDoc.createElement("graphml")
        # define namespaces (proper XML and new expat needs it)
        graphml.setAttribute('xmlns', 'http://graphml.graphdrawing.org/xmlns')
        graphml.setAttribute('xmlns:gd', 'gdnamespace') # find the correct URI
        graphml.setAttribute('xmlns:hmm', 'http://www.ghmm.org/xml/')#arbitrary
        XMLNode.appendChild(graphml)

        # Create key elements
        hmmtype = XMLDoc.createElement("key")
        hmmtype.setAttribute('id', 'emissions')
        hmmtype.setAttribute('gd:type', 'HigherDiscreteProbDist') # what's your type?
        hmmtype.setAttribute('for', 'node')
        graphml.appendChild(hmmtype)
        
        self.hmmClass.toDOM(XMLDoc, graphml)

        if (self.modelType == "pairHMM"):
            modelType = XMLDoc.createElement("hmm:modeltype")
            modelType.appendChild(XMLDoc.createTextNode("pairHMM"))
            graphml.appendChild(modelType)
        
        for alphabet in self.hmmAlphabets.values():
            alphabet.toDOM(XMLDoc, graphml)
        self.backgroundDistributions.toDOM(XMLDoc, graphml) 

        if len(self.transitionFunctions.keys()) != 0:
            transitionFunctionsNode = XMLDoc.createElement("hmm:transitionfunctions")
            for transitionFunction in self.transitionFunctions.values():
                transitionFunction.toDom(XMLDoc, transitionFunctionsNode)
            graphml.appendChild(transitionFunctionsNode)

        graph = XMLDoc.createElement("graph")

        # Compute sums of initial probabilities for renormalization 
        initial_sum = 0.0
        for s in self.state.keys():
            initial_sum = initial_sum + self.state[s].initial
        
        for s in self.state.keys():
            self.state[s].toDOM(XMLDoc, graph, initial_sum)
        
        # Compute sums of outgoing probabilities for renormalization of transition probabilities
        # NOTE: need dictionaries here
        out_sum = {}
        nr_classes = int(self.hmmClass.high())-int(self.hmmClass.low())+1
        for v in self.G.vertices:
            out_sum[v] = [0.0]*nr_classes

        for cl in range(1): # XXX Assuming one transition class
            for e in self.G.Edges():
                if self.G.edgeWeights[cl].has_key(e):
                    out_sum[e[0]][cl] = out_sum[e[0]][cl] + self.G.edgeWeights[cl][e]
                
        for e in self.G.Edges():
            transitions = []
            edge_elem = XMLDoc.createElement("edge")
            edge_elem.setAttribute('source', "%s" % self.state[e[0]].id)
            edge_elem.setAttribute('target', "%s" % self.state[e[1]].id)
            # writeData(XMLDoc, edge_elem, 'prob', self.G.edgeWeights[cl][e] / out_sum[e[0]])
            # XXX Assuming one transition class for cl in range(nr_classes):
            for cl in range(1):
                if self.G.edgeWeights[cl].has_key(e) and out_sum[e[0]][cl]:
                    transitions.append(self.G.edgeWeights[cl][e]/ out_sum[e[0]][cl])
                else:
                    transitions.append(0.0)
                
            writeData(XMLDoc, edge_elem, 'prob', csvFromList( transitions ))

            graph.appendChild(edge_elem)  
            
        graphml.appendChild(graph)

    def AlphabetType(self):
	""" return the type of emission domain 
	    XXX should call the method in HMMAlphabet
	"""
	return int
    
    def ClassType(self):
	pass
    
    def DistributionType(self):
	pass

    def getBackgroundDist(self):
        """ Return a pair of two dictionaries: (distribution, its orders):
            a distribution is a list of real values of length N^(order+1).   
        """        
        return (self.backgroundDistributions.dist, self.backgroundDistributions.order, self.backgroundDistributions.code2name)

        
    def buildMatrices(self):    
	""" return [alphabets_code, A, B, pi, state_orders] """
	pi = []
	B  = []
	A  = []
	nstates = len(self.state.keys())
	orders = {}
	k = 0 # C style index
	for s in self.state.values(): # ordering from XML
	    orders[s.index] = k
	    k = k + 1

        state_orders = []
	for s in self.state.values(): # a list of indices
	    pi.append(s.initial)
	    state_orders.append(s.order) # state order

            size = self.hmmAlphabets[s.alphabet_id].size()
            if (self.modelType != "pairHMM" and
                size**(s.order+1) != len(s.emissions)): 
		raise ValueError # exception: inconsistency between ordering and emission
	    
            B.append(s.emissions) # emission
	    
	    # transition probability
	    v = s.index
	    outprobs = [0.0] * nstates
	    for outid in self.G.OutNeighbors(v)[:]:
		myorder = orders[outid]
		outprobs[myorder] = self.G.edgeWeights[0][(v,outid)]
	    A.append(outprobs)

        alphabets = self.hmmAlphabets[0].name.values() # list of alphabets
	return [alphabets, A, B, pi, state_orders]

    def getStateAlphabets(self):
        alphabets = []
        for s in self.state.values():
            alphabets.append(self.hmmAlphabets[s.alphabet_id])
        return alphabets

    def getAlphabets(self):
        return self.hmmAlphabets
    
    def getLabels(self):
        """ returns list of state labels and unique labels """
        label_list = []
        labels = {}
        for s in self.state.values(): # a list of indices
           label_list.append(self.hmmClass.code2name[s.state_class])
           labels[label_list[-1]] = 0
        return (label_list, labels.keys())

    def getTiedStates(self):    
        """ returns list of tied states, entry is None if a state isn't to
            any other state, returns an empty list, if no state is tied """
        tiedstates = []
        isTied = 0
        
        orders = {}
        k = 0 # C style index
        for s in self.state.values(): # ordering from XML
            orders[s.id] = k
            k = k + 1

        for s in self.state.values(): # a list of indices
            if s.tiedto == '':
                tiedstates.append(-1)
            else:
                tiedstates.append(orders[int(s.tiedto)])
                isTied = 1

        if not isTied:
            tiedstates = []
        return tiedstates

    def getStateDurations(self):
        """ returns a list of the minimal number of times a state is evaluated
            before the HMM changes to another state."""

        durations = []
        hasduration = 0
        
        for s in self.state.values(): # a list of indices
            if s.duration==0:
                durations.append(1)
            else:
                durations.append(s.duration)
                hasduration = 1

        if not hasduration:
            durations = []
        return durations
    
    def OpenXML(self, fileName_file_or_dom):
        if (not isinstance(fileName_file_or_dom, xml.dom.minidom.Document)):
            dom = xml.dom.minidom.parse(fileName_file_or_dom)
        else:
            dom = fileName_file_or_dom
        if dom.documentElement.tagName == "ghmm":
            sys.stderr.write("Do not support ghmm format")
            raise FormatError
            dom.unlink()
            #self.DocumentName = "ghmm"
            #ghmmdom  = dom
            #ghmml = GHMMXML()
            #dom   = ghmml.GraphMLDOM(ghmmdom)
            #ghmmdom.unlink()
        else:
            assert dom.documentElement.tagName == "graphml"   
	    self.fromDOM(dom)
	    # dom.unlink()

    def WriteXML(self, fileName):
        try:
            self.modelCheck()   # raise exceptions here
            doc = xml.dom.minidom.Document()
            self.toDOM(doc, doc)
            file = open(fileName, 'w')
            # xml.dom.ext.PrettyPrint(doc, file)        
            file.write(toprettyxml(doc)) # problem with white spaces
            file.close()
            doc.unlink()
        except HMMEdError:
            print "HMMEdError: No file was written due to errors in the model."
            
    def WriteGHMM(self, fileName):
	self.modelCheck()   # raise exceptions here
        doc = xml.dom.minidom.Document()
        ghmm = doc.createElement("ghmm")
        doc.appendChild(ghmm)
        self.toGHMM(doc, ghmm)
        file = open(fileName, 'w')
        # xml.dom.ext.PrettyPrint(doc, file)        
        file.write(toprettyxml(doc)) # problem with white spaces
        file.close()
        doc.unlink()
        
    def SaveAs(self, fileName):
        if ( self.DocumentName == "graphml" ):
            self.WriteXML(fileName)
        else:
            self.WriteGHMM(fileName)
            
    def SaveAsGHMM(self, fileName):
        self.WriteGHMM(fileName)

class TransitionFunction:
    """ this class holds information on the function which determines the
        transition class for the state and holds the necesary parameters """

    def __init__(self, type=None, paramDict=None):
        self.type = type
        self.paramDict = paramDict

    def fromDom(self, XMLNode):
        self.name = XMLNode.getAttribute("hmm:name")
        self.id = int(XMLNode.getAttribute("hmm:id"))
        typeNode = XMLNode.getElementsByTagName("hmm:transitiontype")[0]
        self.type = typeNode.firstChild.nodeValue
        parameterNodes = XMLNode.getElementsByTagName("hmm:transitionparameter")
        self.paramDict = {}
        for parameter in parameterNodes:
            self.paramDict[parameter.getAttribute("hmm:name")] = parameter.firstChild.nodeValue

    def toDom(self, XMLDoc, XMLNode):
        transitionFunctionNode = XMLDoc.createElement("hmm:transitionfunction")
        transitionFunctionNode.setAttribute("hmm:name", self.name)
        transitionFunctionNode.setAttribute("hmm:id", self.id)
        XMLNode.appendChild(transitionFunctionNode)
        typeNode = XMLDoc.createElement("hmm:transitiontype")
        type = XMLDoc.createTextNode(self.type)
        typeNode.appendChild(type)
        transitionFunctionNode.appendChild(typeNode)
        for paramName in self.paramDict.keys():
            paramNode = XMLDoc.createElement("hmm:transitionparameter")
            paramNode.setAttribute("hmm:name", paramName)
            param = XMLDoc.createTextNode(self.paramDict[paramName])
            paramNode.appendChild(param)
            transitionFunctionNode.appendChild(paramNode)
        
################################################################################
if __name__ == '__main__':

    hmmobj = HMM()
    hmmobj.OpenXML(sys.argv[1])
    hmmobj.WriteXML("utz.xml")
    

