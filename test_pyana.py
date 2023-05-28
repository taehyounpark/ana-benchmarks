import ROOT

from ROOT import ana
import cppyy

ds = ROOT.ana.analysis['Tree']([],'')
# one = ds.constant(1)
# cppyy.cppdef(''' auto plus_one_fn = std::function([](double x){return x+1;}); ''')
# plus_one = ds.define(cppyy.gbl.plus_one_fn)
# two = plus_one.evaluate(one)  # <- this will result in error

# print(ds)
# print(one)
# print(cppyy.gbl.plus_one_fn)
# print(plus_one)

# directly to cling
cppyy.cppdef(''' auto ds = ana::analysis<Tree>({},""); ''')
print(cppyy.gbl.ds)
cppyy.cppdef(''' auto one = ds.constant(1); ''')
print(cppyy.gbl.one)
cppyy.cppdef(''' auto plus_one = ds.define([](double x){return x+1;}); ''')
print(cppyy.gbl.plus_one)
ds = cppyy.gbl.ds
one = cppyy.gbl.one
plus_one = cppyy.gbl.plus_one
ds.evaluate_column(plus_one,one)
# plus_one.evaluate_column(one)
# cppyy.cppdef(''' auto two = plus_one(one); ''')

# print(cppyy.gbl.two)
