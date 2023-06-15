'''
File di appoggio per il corso di Python

Sections:

> Mutable e immutable       -> dummy_muta
> int                       -> dummy_int
> Factorial                 -> dummy_fact


'''

import numpy as np

def dummy_muta():
    '''
    '''
    print("Esempio immutable:")
    a=5
    b=a
    print('a = ',a,' b =',b)
    print('id(a)',id(a),'id(b)',id(b),'id(5)',id(5))
    a=34
    print('a = ',a,' b =',b)
    print('id(a)',id(a),'id(b)',id(b),'id(5)',id(5),'id(34)',id(34))

    print("\nEsempio mutable:")
    l=[1,2,3]
    l1=l
    print(l,l1)
    print('id(l)',id(l),'id(l1)',id(l1))
    l[1]=34
    print(l,l1)
    print('id(l)',id(l),'id(l1)',id(l1))



def dummy_int():
    '''
    '''
    a=1000000
    a=a*1000000
    a=a*1000000
    a=a*1000000
    a=a*1000000
    a=a*1000000
    a=a*1000000
    a=a*1000000
    a=a*1000000
    a=a*1000000
    a=a*1000000
    a=a*1000000
    a=a*1000000
    a=a*1000000
    a=a*1000000
    a=a*1000000
    a=a*1000000

    print(a)
    print(float(a))
    b=a+345
    print('a,b,b-a',a,b,b-a)



def dummy_fact(n):
    '''
    '''
    print('   Uso di ifact')
    ifact(n)

    print('\n   Uso di lffact')
    lffact(n)

    print(' \n  Uso di rec_fact')
    r=rec_fact(n)
    print('The factorial of',n,'is',r)

    print('\n   Uso di Stirling')
    stirling(n)

    print('\n   Uso di ffact')
    ffact(n)



def ifact(n):
    '''
    semplice
    '''
    fact = 1
    
    for i in range(1, n+1):
        fact = fact * i
    
    print('The factorial of',n,' is : ', end="")
    print(fact) 

    
def ffact(n):
    '''
    semplice in float
    '''
    fact = 1.
    
    for i in range(1, n+1):
        fact = fact * i
    
    print('The factorial of',n,' is : ', end="")
    print(fact) 


def rec_fact(n):
    '''
    semplice ricorsiva
    '''
    if n < 2:
        return 1
    else:
        return n * rec_fact(n-1)
    

def lffact(n):
    '''
    log float
    '''

    l10=np.log10(np.arange(1,n+1))

    A=sum(l10)
    IA=int(A)
    a=10**(A-IA)

    print(a,'10^',IA)


def stirling(n):
    '''
    formula di Stirling
    '''
    
    r=np.sqrt(2*np.pi*n)*(n/np.e)**n

    print('Formula di Stirling',r)

