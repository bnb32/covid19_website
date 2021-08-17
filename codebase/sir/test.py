args={'arg':20}

def test2(**args):
    args['t']=30
    print(test(**args))

def test(t,arg=10):
    print(t,arg)
    return 0

print(test2(**args))    
