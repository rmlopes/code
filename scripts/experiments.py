rates = {'_14': '0.1,0.4',
         '_23': '0.2,0.3',
         '_32': '0.3,0.2',
         '_41': '0.4,0.1',
         '_55': '0.25,0.25'}

blind_operators = {'TD': 'transposon,delete',
                   'JD': 'junk,delete'}

sizes = {'_50': '50',
         '_100': '100',
         '_200': '200',
         '_400': '400'}

aware_operators = 'genecopy,genedelete'

xoverops = {'1P' : 'onepoint',
            '2P' : 'twopoint',
            'uni': 'uniform',
            '1Pgene': 'onepointgene',
            '2Pgene': 'twopointgene',
            'unigene': 'unigene'}

xo_rates = ['0.1','0.3','0.5','0.7','0.9']

#Variation operators experiments
ops = dict()

for oplabel, opv in blind_operators.items():
    for slabel,s in sizes.items():
        for rlabel, rv in rates.items():
            l = [('--operators '+ opv)]
            l.append('--opargs ' + s)
            l.append('--oprates ' + rv)
            ops[oplabel+slabel+rlabel] = l

for rlabel, rv in rates.items():
    l = ['--operators '+ aware_operators]
    l.append('--oprates ' + rv)
    ops['GCD'+rlabel] = l


#Xover operators experiments
xover = dict()
for oplabel, opv in xoverops.items():
    for r in xo_rates:
        l = [('--xover '+ opv + '_xover')]
        l.append('--xrate ' + r)
        xover[oplabel+'-'+r] = l


representation = {'dm5T':['--initdm 5',
                          '--agent code.rencode.DMAgent',
                          '-o True'],
                  'rnd5T':['--initdm 5',
                           '--agent code.rencode.RndAgent',
                           '-o True'],
                  'dm6T':['--initdm 6',
                          '--agent code.rencode.DMAgent',
                          '-o True'],
                  'rnd6T':['--initdm 6',
                           '--agent code.rencode.RndAgent',
                           '-o True'],
                  'dm7T':['--initdm 7',
                          '--agent code.rencode.DMAgent',
                          '-o True'],
                  'rnd7T':['--initdm 7',
                           '--agent code.rencode.RndAgent',
                           '-o True'],
                  'dm8T':['--initdm 8',
                          '--agent code.rencode.DMAgent',
                          '-o True'],
                  'rnd8T':['--initdm 8',
                           '--agent code.rencode.RndAgent',
                           '-o True'],
                  'dm9T':['--initdm 9',
                          '--agent code.rencode.DMAgent',
                          '-o True'],
                  'rnd9T':['--initdm 9',
                           '--agent code.rencode.RndAgent',
                           '-o True'],
                  'dm6F':['--initdm 6',
                          '--agent code.rencode.DMAgent',
                          '-o False'],
                  'rnd6F':['--initdm 6',
                           '--agent code.rencode.RndAgent',
                           '-o False'],
                  'dm7F':['--initdm 7',
                          '--agent code.rencode.DMAgent',
                          '-o False'],
                  'rnd7F':['--initdm 7',
                           '--agent code.rencode.RndAgent',
                           '-o False'],
                  'dm8F':['--initdm 8',
                          '--agent code.rencode.DMAgent',
                          '-o False'],
                  'rnd8F':['--initdm 8',
                           '--agent code.rencode.RndAgent',
                           '-o False'],
                  'dm9F':['--initdm 9',
                          '--agent code.rencode.DMAgent',
                          '-o False'],
                  'rnd9F':['--initdm 9',
                           '--agent code.rencode.RndAgent',
                           '-o False']}

if __name__ == '__main__':
    print ops
