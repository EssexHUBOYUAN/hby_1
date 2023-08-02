# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 11:09:32 2023

this project aims to generate test sample blocks, and corresponding ingredients given corrilation between rate and ing across age.

@author: Doc Who
"""
import csv
from scipy.optimize import curve_fit
import numpy as np
from scipy import special as sp
import scipy as scp
from scipy import stats
import math
import matplotlib.pyplot as plt
import random

'''form population file, given total population of a year'''
def get_population_age(pop, oldage, deathrange):
    arr_age = np.arange(120)
    arr_pop = np.array([StepEXP(age,0-deathrange,oldage) for age in arr_age])
    arr_pop = pop/np.sum(arr_pop)*arr_pop
    arr_pop_86 = list(arr_pop[:85]) + [np.sum(arr_pop[85:])]
    return [round(a) for a in arr_pop_86]

def StepEXP(x,a,b):
    '''
    if a > 0:
        print('StepEXP is 0 --> 1')
    else:
        print('StepEXP is 1 --> 0')
    '''
    return np.exp((x-b)/a)/(1+np.exp((x-b)/a))

'''___________________________________________________________________________________________________________________'''
'''get the parameters for expected average rate at given age (0-85+)'''

L_initis = [[0.05, 0.0001, 66, 7],[0.01, 0.0001, 66, 7],[0.005, 0.00006, 66, 7],[0.005, 0.00006, 63, 7]] # for SEER17
def read_summary_rate_file(filename):
    D_data = dict()
    with open(filename+'.csv', mode='r',encoding='UTF-8-sig') as file:
        data = csv.reader(file)
        Head = True
        for line in data:
            if Head:
                Head = False
            else:
                D_data[line[0]] = [float(line[i+1]) for i in range(len(line)-1)]
    D_params = dict()
    lx = [5*i for i in range(18)]
    for gender in D_data:
        ly = D_data[gender]
        D_params[gender] = age_fit(lx,ly)
    return D_params
def age_curve(x, a,b,c,d):
    return -a/b*(np.log(np.sqrt(b)*x+1)-np.sqrt(b)*x)/(np.exp((x-c)/d)+1)
def age_fit(lx,ly):
    for candidate in L_initis:
        TryNext = False
        NotBest = True
        try:
            initis = candidate
            fitParams, fitCovariances = curve_fit(age_curve, lx, ly, initis)
            a,b,c,d = fitParams
            if a < 0  or c > 125 or c < 12 or d < 0:
                #print('try next initis for age curve fitting')
                TryNext = True
            else:
                pass
            if TryNext:
                pass
            else:
                #print(fitParams,' is found')
                NotBest = False
                break
        except RuntimeError:
            #print('try next initis for age curve fitting')
            pass
    if NotBest:
        print('finally', fitParams,' is found')
    else:
        pass
    #print(fitParams)
    return fitParams
'''__________________________________________________________________________________________________________________'''
'''get a list of corr_random, for a given value of corr_exp'''
def get_corrs(corr, N_countrys, N_years):
    lx = [i*0.01-1 for i in range(200)]
    ly0 = [corr_distribution(x, N_years, corr) for x in lx]
    L_choosed = random.choices(lx,ly0,k = N_countrys)
    L_choosed.sort()
    return L_choosed

def SUM(List):
    S = 0
    for item in List:
        S = S + item
    return S

def corr_distribution(x, n, lo=0):
    if lo == 0:
        return (1-x**2)**((n-1)/2)*(sp.gamma(n/2)/sp.gamma((n-1)/2)/sp.gamma(1/2))
    else:
        L = [(2*lo*x)**k/math.factorial(k)*(sp.gamma((n+k)/2))**2 for k in range(35)]
        return abs(2**(n-2)*(1-lo**2)**(n/2)*(1-x**2)**((n-1)/2-1)/(math.factorial(n-2)*math.pi)*SUM(L))

def cum_corr_distribution(x,n,lo=0):
    lx = []
    a = -1
    while a < x:
        a = a+0.01
        lx.append(a)
    ly = [corr_distribution(a, n, lo)/100 for a in lx]
    return SUM(ly)
'''__________________________________________________________________________________________________________________'''
'''generate annual average rate vs ingredient. 2 lists of N_years'''
def get_rates_vs_ing(ybar, dy, xbar, dx, corr, N_years):
    if corr == 0 or dx == 0 or dy == 0:
        slope = 0
    else:
        slope = ( (dx**2-dy**2) - np.sqrt( (dx**2-dy**2)**2 + 4*(corr*dy*dx)**2 ) )/(-2*corr*dy*dx)
    da = np.sqrt( (1+slope**2)**(-1)*(dx**2 + slope**2*dy**2 + 2*slope*corr*dx*dy) )
    db = np.sqrt( (1+slope**2)**(-1)*(dy**2 + slope**2*dx**2 - 2*slope*corr*dx*dy) )
    abar = (xbar + slope*ybar)/np.sqrt(1+slope**2)
    bbar = (-slope*xbar + ybar)/np.sqrt(1+slope**2)
    l_a = np.random.normal(abar, da, N_years)
    l_b = np.random.normal(bbar, db, N_years)
    arr_pairs = np.transpose(np.array([l_a,l_b]))
    l_x = []
    l_y = []
    for pair in arr_pairs:
        a, b = pair
        delta = a-abar
        deltb = b-bbar
        l_x.append( (delta - deltb*slope)/(1+slope**2)+xbar )
        l_y.append( (delta*slope + deltb)/(1+slope**2)+ybar )
    '''
    plt.title('this fig has corr: '+str(round(corr,2)) )
    plt.scatter(l_x, l_y)
    plt.scatter(SUM(l_x)/len(l_x), SUM(l_y)/len(l_y))
    plt.scatter(xbar,ybar)
    plt.show()
    '''
    return np.transpose(np.array([l_x,l_y]))

'''__________________________________________________________________________________________________________________'''
'''form counts block shape:(86,Nyears), then shape:(18,N_years), and rate-raw:(18,N_years), pop:(18,N_years)'''
def count_random_possion(p_exp, N):
    la = N*p_exp/10**5
    if la < 24:
        lweight = [stats.poisson.pmf(x,la) for x in range(int(10*la)+1)]
        return random.choices([x for x in range(int(10*la)+1)], lweight, k=1)[0]
    else:
        std = np.sqrt(la)
        return random.normalvariate(la, std)
def Possion_dis(la, x):
    if int(x) == float(x):
        k = int(x)
        Ga = 1
        for i in range(k):
            Ga = Ga*(i+1)
            if Ga > 7257415615307998967396728211129263114716991681296451376543577798900561843401706157852350749242617459511490991237838520776666022565442753025328900773207510902400430280058295603966612599658257104398558294257568966313439612262571094946806711205568880457193340212661452800000000000000000000000000000000000000000:
                Ga = 7257415615307998967396728211129263114716991681296451376543577798900561843401706157852350749242617459511490991237838520776666022565442753025328900773207510902400430280058295603966612599658257104398558294257568966313439612262571094946806711205568880457193340212661452800000000000000000000000000000000000000000
            else:
                pass
        p = np.exp(-la)*la**k/Ga
        #print('Ga:',Ga)
    else:
        p = "ValueError: x not int"
    
    return p

def get_counts86( l_pops, arr2_rate86):
    a, b = np.shape(arr2_rate86)
    # a = number of years, b = age 0~85
    arr2_counts = np.zeros((a,b))
    for i in range(a):
        for j in range(b):
            arr2_counts[i,j] = count_random_possion(arr2_rate86[i,j], l_pops[j])
    return arr2_counts

def get_age5(l_pops):
    L = []
    for i in range(17):
        l = []
        for j in range(5):
            ind = 5*i + j
            l.append(l_pops[ind])
        a = SUM(l)
        L.append(a)
    L.append(l_pops[-1])
    return L

def get_counts5(arr2_counts):
    L = []
    for year_counts86 in arr2_counts:
        year_counts = get_age5(year_counts86)
        L.append(year_counts)
    return np.array(L)

def get_pops5(N_years, l_pops86):
    L = []
    for i in range(N_years):
        l_pops5 = get_age5(l_pops86)
        L.append(l_pops5)
    return np.array(L)

def get_rate5(arr_counts, arr_pops):
    return arr_counts/arr_pops*10**5

'''__________________________________________________________________________________________________________________'''
'''save different files'''
def save_pops_file(D_pop, country, L_years):
    #print(D_pop['male'])
    a,b = np.shape(D_pop['male'])
    if b == 86:
        headline = ['region','gender','year'] + [i for i in range(b)]
        filename = country + '_pop_86org.csv'
    else:
        headline = ['region','gender','year'] + [i*5 for i in range(b)]
        filename = country + '_pop.csv'
    with open(filename, mode='w', newline = '') as file:
        data = csv.writer(file)
        data.writerow(headline)
        for gender in D_pop:
            for i in range(a):
                L = []
                for pop in D_pop[gender][i]:
                    L.append(pop)
                data.writerow([country, gender, L_years[i]] + L)

    return print(country, '_pop file saved, with age length: ', b)

def save_counts_file(D_counts, country, L_years):
    a,b = np.shape(D_counts[country,'male'])
    if b == 86:
        headline = ['region','gender','year'] + [i for i in range(b)]
        filename = country + '_counts_86org.csv'
    else:
        headline = ['region','gender','year'] + [i*5 for i in range(b)]
        filename = country + '_counts.csv'
    with open(filename, mode='w', newline = '') as file:
        data = csv.writer(file)
        data.writerow(headline)
        for item in D_counts:
            if item[0] == country:
                for i in range(a):
                    L = []
                    for count in D_counts[country, item[1]][i]:
                        L.append(count)
                    data.writerow([country, item[1], L_years[i]] + L)
            else:
                pass
    return print(country, '_counts file saved, with age length: ', b)

def save_rate_file(D_rate, country, L_years):
    a,b = np.shape(D_rate[country,'male'])
    if b == 86:
        l_age = [i for i in range(b)]
        filename = country + '_rate_org.csv'
    else:
        l_age = [i*5 for i in range(b)]
        filename = country + '_raw.csv'
    headline1 = ['','']
    headline2 = ['blocksize', 'age\year']
    L_blocks = []
    for item in D_rate:
        gender = item[1]
        if item[0] == country:
            headline1 = headline1 + [gender for i in range(a)]
            headline2 = headline2 + L_years
            L_blocks.append(D_rate[country,gender])
        else:
            pass
    block = np.concatenate([np.transpose(L_blocks[0]), np.transpose(L_blocks[1])], axis=1)
    with open(filename, mode='w', newline = '') as file:
        data = csv.writer(file)
        data.writerow(headline1)
        data.writerow(headline2)
        FIRSTLINE = True
        for row in range(b):
            if FIRSTLINE:
                line = ['0X0',0] + [rate for rate in block[row]]
                FIRSTLINE = False
            else:
                line = ['',l_age[row]] + list(block[row])
            data.writerow(line)
    return print(country, '_rates file saved, with age length: ', b)

def save_ingredients(D_ing, country, N_years):
    filename = country + '_ingredients.csv'
    headline = ['Geography','Category'] + N_years
    with open(filename, mode='w', newline='') as file:
        data = csv.writer(file)
        data.writerow(headline)
        for item in D_ing:
            if item[0] == country:
                product = item[1]
                line = [country, product] + [10**x for x in D_ing[country,product]]
                data.writerow(line)
            else:
                pass
    return print(country, '_ingredients file saved' )

'''__________________________________________________________________________________________________________________'''
'''main'''

def generate_files(L_corr_exp, N_sims, filehead):
    D_ingredients = dict() #country, product: list of ing(in log-tonne)
    D_rate5 = dict() # male/female: rate block (age_ave5 vs years)
    D_pops5 = dict() # male/female: pop block (age_sum5 vs years)
    D_counts5 = dict() #male/female: counts block (counts_sum5 vs years)
    D_rate86 = dict() # male/female: rate block (age_ave86 vs years)
    D_pops86 = dict() # male/female: pop block (age_86 vs years)
    D_counts86 = dict() #male/female: counts block (counts_sum86 vs years)
    L_country = set()
    #L_products = set()
    
    for gender in D_params:
        a,b,c,d = D_params[gender]
        print(gender, "old ages:", round(c,2),round(d,2))
        ly = [age_curve(age,a,b,c,d) for age in lx]
        pops = get_population_age(N_sims, c, d)
        pops5 = get_pops5(N_years, pops)
        D_pops86[gender] = np.array([pops for i in range(N_years)])
        print('D_pops86',np.shape(D_pops86[gender]))
        D_pops5[gender] = np.array(pops5)
        print('D_pops5',np.shape(D_pops5[gender]))
        ave_rate = np.sum(np.array(ly)*(np.array(pops)/np.sum(np.array(pops))))
        print(gender, ave_rate)
        #plt.plot(lx,pops)
        product_ind = 1
        for corr_exp in L_corr_exp:
            product = filehead +'product'+ str(product_ind) + '&' + gender
            #L_products.add(product)
            product_ind = product_ind + 1
            L_corrs = get_corrs(corr_exp, N_countrys, N_years)
            print('number of countrys: ', len(L_corrs))
            arr2_pops5 = get_pops5(N_years, pops)
            country_ind = 1
            for corr in L_corrs:
                country = filehead + 'country' + str(country_ind)
                L_country.add(country)
                country_ind = country_ind + 1
                # each corr refer to a country
                arr_pairs_eachyear = get_rates_vs_ing(ave_rate, std_rate_frac*ave_rate, ave_ing, std_ing, corr, N_years)
                arr2_rate_age_year = []
                arr1_ing_years = []
                for pairs in arr_pairs_eachyear:
                    # each pair refer to a year
                    rate_year = pairs[1]
                    arr1_ing_years.append(pairs[0])
                    arr2_rate_age_year.append( np.array(ly)/np.sum(np.array(ly))*rate_year*86 )
                
                arr2_rate_age_year = np.array(arr2_rate_age_year)
                D_rate86[country, gender] = arr2_rate_age_year
                D_ingredients[country, product] = arr1_ing_years
                
                arr2_counts86 = get_counts86(pops, arr2_rate_age_year)
                D_counts86[country, gender] = arr2_counts86
                
                arr2_counts5 = get_counts5(arr2_counts86)
                D_counts5[country, gender] = arr2_counts5
                
                arr2_rate_raw = get_rate5(arr2_counts5, arr2_pops5)
                D_rate5[country, gender] = arr2_rate_raw
    L_country = list(L_country)
    #L_products = list(L_products)
    # save pops, counts, and raw_rate file (18 and 86)   
    for country in L_country:
        save_counts_file(D_counts86, country, L_years)
        save_counts_file(D_counts5, country, L_years)
        save_pops_file(D_pops86, country, L_years)
        save_pops_file(D_pops5, country, L_years)
        save_rate_file(D_rate86, country, L_years)
        save_rate_file(D_rate5, country, L_years)    
        # save ing file for different ings 
        save_ingredients(D_ingredients, country, L_years)   
    #plt.show()

lx = [age for age in range(86)]
std_ing = 1
ave_ing = 3
std_rate_frac = 0.3
N_years = 24
L_years = [(i+2000) for i in range(N_years)]
N_countrys = 24
filename = 'country_summary_rates'
D_params = read_summary_rate_file(filename)

def test_generator(list_corrs = [[0.1],[0.2],[0.3],[0.4],[0.5],[0.6],[0.7],[0.8],[0.9]], list_N_sims = [5*10**5, 10**6, 5*10**6, 10**7, 5*10**7, 10**8, 5*10**8, 10**9], model_index = 0):
    for corr_exp in list_corrs:
        corr_exp = [corr_exp[0] + 0.0000000001*model_index]
        for N_sims in list_N_sims:
            filehead = 'cor'+str(round(corr_exp[0],1))+'^'+'pop'+str(round(np.log10(N_sims),1))+'^'
            generate_files(corr_exp, N_sims, filehead)
    return