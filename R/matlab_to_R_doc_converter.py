# Convert MATLAB documentation to R Rd doc files

import re, sys, subprocess, os
from optparse import OptionParser

keywords = ['DESC', 'ARG', 'FORMAT', 'RETURN', 'SEEALSO', 'COPYRIGHT', 'MODIFICATIONS']

BASEPATH = os.path.expanduser('~/mlprojects/.')
MANPATH = os.path.expanduser('~/mlprojects/./gpsim/R/man')

R_FUNC_RE = re.compile(r"(\w+) <- function")


def find_matlab_file(basepath, func):
    p = subprocess.Popen(['find',basepath,'-name','%s.m' % func], stdout=subprocess.PIPE)
    p.wait()
    l = p.stdout.next().strip()
    return l


def find_R_functions(files):
    v = []
    for fn in files:
        f = open(fn)
        for l in f:
            m = R_FUNC_RE.match(l)
            if m:
                v.append(m.group(1))
        f.close()
    return v


def matlab_doc_continuation(f, text):
    while True:
        l = f.next()
        if l[0:2] != "% ":
            return (l, text)
        t = l[2:].partition(' ')
        if t[0] not in keywords:
            text = text + ' ' + l[2:].strip()
        else:
            return (l, text)
    


def parse_matlab_doc(fname, func):
    f = open(fname)
    doc = {'name': func, 'formats': []}
    #'args': [], 'returns': [], 'desc': []}
    l = f.next()
    while l[0] != "%":
        l = f.next()
    
    l = l[2:].strip()
    if l.startswith(func.upper()):
        t = l.partition(' ')
        doc['title'] = t[2]
    l = f.next()

    while True:
        #print l
        if l[0] != "%":
            break
        l = l[2:].strip()
        if l.startswith("FORMAT"):
            doc['formats'].append({'args': [], 'returns': []})
            #doc['args'].append([])
            #doc['returns'].append([])
            #doc['desc'].append([])
            l = f.next()
        elif l.startswith("DESC"):
            t = l.partition(' ')
            (l, desc) = matlab_doc_continuation(f, t[2])
            doc['formats'][-1]['desc'] = desc
        elif l.startswith("ARG"):
            t = l[3:].partition(':')
            (l, desc) = matlab_doc_continuation(f, t[2].strip())
            doc['formats'][-1]['args'].append([t[0].strip(), desc])
        elif l.startswith("RETURN"):
            t = l[6:].partition(':')
            (l, desc) = matlab_doc_continuation(f, t[2].strip())
            doc['formats'][-1]['returns'].append([t[0].strip(), desc])
        elif l.startswith("SEEALSO"):
            t = l.partition(':')
            (l, desc) = matlab_doc_continuation(f, t[2].strip())
            doc['seealso'] = desc.split(',')
        elif l.startswith("COPYRIGHT"):
            t = l.partition(':')
            (l, desc) = matlab_doc_continuation(f, t[2].strip())
            doc['copyright'] = desc
        elif l.startswith("MODIFICATIONS"):
            t = l.partition(':')
            (l, desc) = matlab_doc_continuation(f, t[2].strip())
            doc['modifications'] = desc
        else:
            l = f.next()
    f.close()
    return doc


def linewrap_R_call(s):
    # split arguments
    t = s.split(',')
    # length if the initial part, use to indent next lines
    clen = s.index('(')
    # lines to return
    v = [t[0]]
    for k in t[1:]:
        # next token fits this line:
        if len(','.join([v[-1], k])) < 65:
            v[-1] = ','.join([v[-1], k])
        # next token needs new line:
        else:
            v[-1] += ','
            v.append(' '*clen + k)
    return '\n'.join(v)


def make_R_call(name, format, valuelist=False):
    # "err = modelObjective(params, model, \dots)"
    arglist = ', '.join(map(lambda x: x[0], format['args']))
    if len(format['returns']) == 0:
        val = "%s(%s)" % (name, arglist)
    elif len(format['returns']) == 1:
        val = "%s <- %s(%s)" % (format['returns'][0][0], name, arglist)
    else:
        if valuelist:
            vallist = ', '.join(map(lambda x: x[0], format['returns']))
            val = "%s <- %s(%s)" % (vallist, name, arglist)
        else:
            val = "values <- %s(%s)" % (name, arglist)
    if len(val) > 65:
        return linewrap_R_call(val)
    else:
        return val


def write_R_doc(doc, fout=sys.stdout):
    print >>fout, "\\name{%s}" % doc['name']
    print >>fout, "\\Rdversion{1.0}"
    print >>fout, "\\alias{%s}" % doc['name']
    print >>fout, "\\title{%s}" % doc['title']
    print >>fout, "\\description{"
    print >>fout, "  %s" % doc['formats'][0]['desc']
    print >>fout, """}
\usage{"""
    print >>fout, make_R_call(doc['name'], doc['formats'][0])
    print >>fout, """}
\\arguments{"""
    for k in doc['formats'][0]['args']:
        print >>fout, "  \\item{%s}{%s}" % tuple(k)
    print >>fout, """}
\\value{"""
    for k in doc['formats'][0]['returns']:
        print >>fout, "  \\item{%s}{%s}" % tuple(k)
    print >>fout, """}
\\seealso{"""
    print >>fout, "\\code{" + ', '.join(map(lambda x: '\\link{%s}' % x.strip(), doc['seealso'])) + "}."
    print >>fout, """}
\\examples{
## missing
}
\\keyword{model}"""


def make_args_unique(formats):
    seen_args = dict()
    for i in range(len(formats)):
        for j in range(len(formats[i]['args'])):
            k = formats[i]['args'][j]
            if k[0] in seen_args:
                if seen_args[k[0]] != k[1]:
                    arg = k[0]
                    value = k[1]
                    while arg in seen_args and (seen_args[arg] != value):
                        arg += '\''
                    formats[i]['args'][j][0] = arg
            else:
                seen_args[k[0]] = k[1]
    return formats



def write_R_doc_multiformat(doc, fout=sys.stdout):
    print >>fout, "\\name{%s}" % doc['name']
    print >>fout, "\\Rdversion{1.0}"
    print >>fout, "\\alias{%s}" % doc['name']
    print >>fout, "\\title{%s}" % doc['title']
    print >>fout, "\\description{"
    print >>fout, "  %s" % doc['title']
    print >>fout, """}
\usage{"""
    for f in doc['formats']:
        print >>fout, make_R_call(doc['name'], f)
    print >>fout, """}
\\arguments{"""
    done_args = dict()
    for f in doc['formats']:
        for k in f['args']:
            if not k[0] in done_args:
                print >>fout, "  \\item{%s}{%s}" % tuple(k)
                done_args[k[0]] = k[1]
    print >>fout, """}
\\details{"""
    for f in doc['formats']:
        print >>fout, "  \\code{" + make_R_call(doc['name'], f, valuelist=True) + "}"
        print >>fout, " ", f['desc']
        print >>fout, ""
    print >>fout, """}
\\value{"""
    done_values = set()
    for f in doc['formats']:
        for k in f['returns']:
            if not k[0] in done_values:
                print >>fout, "  \\item{%s}{%s}" % tuple(k)
                done_values.add(k[0])
    print >>fout, """}
\\seealso{"""
    print >>fout, "\\code{" + ', '.join(map(lambda x: '\\link{%s}' % x.strip(), doc['seealso'])) + "}."
    print >>fout, """}
\\examples{
## missing
}
\\keyword{model}"""


usage="%prog [options]\n\n" + \
"""Converts Matlab documentation to R Rd files."""

parser = OptionParser(usage=usage, version="%prog 1.0")
parser.add_option("-R", "--fileinput", action="store_true", dest="filein",
                  help="Arguments are R source files, process all functions defined in them.", default=False)
parser.add_option("-o", "--fileoutput", action="store_true", dest="fileout",
                  help="Direct output to files.", default=False)

(options, args) = parser.parse_args()

if options.filein:
    lst = find_R_functions(args)
else:
    lst = args
for func in lst:
    try:
        fname = find_matlab_file(BASEPATH, func)
    except:
        print "No Matlab source found for %s" % func
        continue
        #sys.exit(1)
    doc = parse_matlab_doc(fname, func)
    doc['formats'] = make_args_unique(doc['formats'])
    #print doc
    if len(doc['formats']) > 0:
        if options.fileout:
            print "Writing documentation for %s to file" % doc['name']
            f = open('%s/%s.Rd' % (MANPATH, doc['name']), 'w')
        else:
            f = sys.stdout
        if len(doc['formats']) < 2:
            write_R_doc(doc, f)
        else:
            write_R_doc_multiformat(doc, f)
    else:
        print "No documentation for %s" % func
