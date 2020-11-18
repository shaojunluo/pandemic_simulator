from numpy import random as rd
from collections import Counter
from scipy.stats import poisson
import networkx as nx
import pandas as pd
import numpy as np

default_params = {
    # Fixed disease parameter
    'b' :  0.0155, # probability of in-community infected rate by symptomatic individuals (tuned)
    'ej' : 0.36,  # factor of getting infected from isolated individuals (strength of quarantine) from [1]
    'p0' : 1e-4,  # probability of getting infected from outside environment
    'DI' : 5,     # median  days of incubation (from WHO)
    'DR' : 10,    # suggest release date after isolation (from WHO)
    'DQ' : 14,    # min days require for quarantine (from NYS)
    # mode parameter
    'lt' : 1.0,   # lockdown threshold
    'ls' : 0.0,    # lockdown scale
    # tunable parameter
    'g'  : 0.0,   # probability of going isolation if shows disease sympton
    'pc' : 0.9,   # protection effectiveness
    'pr' : 0.0,   # rate of person have protection
    'qt' : 0.0,   # the probability of exposion awareness (test rate)
}
    
class person:
    def __init__(self, pid, source = -1, kcore=0, status ='S', infected = False, protected = False,
                       params = default_params):
        self.pid = pid                     # id of person
        self.source = source               # infection source of person if exist
        self.kcore = kcore                 # kcore number of person 
        self.status = status               # status of person 
        self.protected = protected         # if protected measurement is taken
        
        self.disease_status_days = 0       # days of being maintaining disease status
        self.quarantined_days = 0          # days from quarantined
        self.spread = 0                    # how many people it has spred the disease

        self.infected = infected           # Flag showing if you are infeced
        self.lockdown = False              # lockdown status, can only controled by set_lockdown function
        
        # conflict check
        if self.status in ['E','J','I']:
            self.infected = True
        elif self.status in ['S','R']:
            self.infected = False
        elif self.status != 'Q':
            raise ValueError
        
        # assign parameters
        self.b = params['b']
        self.ej = params['ej']
        self.p0 = params['p0']
        self.g = params['g']
        self.DI = params['DI']
        self.DR = params['DR']
        self.DQ = params['DQ']
        self.pc = params['pc']
        self.qt = params['qt']
    
    # count number of status
    def count_status(self,close_contact):
        return Counter([n.status for n in close_contact])
    
    # check if get infected
    def get_infected(self,close_contact):
        counter = self.count_status(close_contact)
        # probability of getting infected
        p = self.p0 + self.b* (self.ej * counter['J'] + counter['I'])
        # add protection
        p = (1-self.pc*int(self.protected))*p
        # infected
        infected = rd.rand() < p
        # assgin infection source
        infected_prob = {'J':self.b*self.ej,'I': self.b}
        if infected:
            candidate = [self] + [c for c in close_contact if c.status in ['J','I']]
            probs = np.array([self.p0] + [infected_prob[c.status] for c in close_contact if c.status in ['J','I']])
            # pick a source and mark the spread
            source = rd.choice(candidate, p=probs/probs.sum())
            source.spread += 1
            source_id = source.pid
        else:
            source_id = -1
        return infected, source_id
    
    # determine wether to go quarantined
    def go_quarantine(self, close_contact,lockdown=False):
        counter = self.count_status(close_contact)
        # number of close contact being confirmed infected
        num_showed_infected = counter['J']+counter['I']
        # if you know you are exposed to a possible sources decide go quarantine
        return self.qt*num_showed_infected > 1
        
    def update(self, close_contact):
        # =================== People with full freedom =================== #
        # suscepitible individuals
        if self.status == 'S':
            # calculate if get infected
            infected, source = self.get_infected(close_contact)
            if infected:
                self.status = 'E'
                self.disease_status_days = 0
                self.infected = True
                self.source = source
                
            # Now check if you need to quaratined
            if self.go_quarantine(close_contact):
                self.status = 'Q'    # go quarantine
                self.quarantined_days = 1
            # if it during lock down, EXTEND the quarantine day!
            if self.lockdown:
                self.status = 'Q'
                self.quarantined_days = -self.DQ
            
        # asymptomatic individuals
        elif self.status == 'E':
            # Chances to developed sympton
            if rd.rand() < poisson.cdf(self.disease_status_days, self.DI):
                self.status = 'I'
                self.disease_status_days = 1
            else:
                self.disease_status_days += 1
                
            # Now check if you need to quaratined
            if self.go_quarantine(close_contact):
                self.status = 'Q'    # go quarantine
                self.quarantined_days = 0 # start quarantine
            # if it during lock down, EXTEND the quarantine day!
            if self.lockdown:
                self.status = 'Q'
                self.quarantined_days = -self.DQ
        
        # if infected but not go to isolated
        elif self.status == 'I':
            # check if go for medical help
            if rd.rand() < self.g:
                self.status = 'J'  # go isolation
                self.disease_status_days += 1 # do not reset status days because we need to count the day of infected
            else:
                self.disease_status_days += 1
                
            # check if recovered
            if self.disease_status_days > self.DR:
                self.status = 'R' # being recovered 
                self.infected = False # mark as no infected
                self.disease_status_days = 0 # reset status days
        
        # if recovered, then nothing to fear nothing to do
        elif self.status == 'R':
            pass
            
        # =================== People with limited freedom =================== #
            
        # if quarantined
        elif self.status == 'Q':
            # add another day of quarantine
            self.quarantined_days += 1
            
            # if you are asymptotic
            if self.infected: 
                # if developed to disease
                if rd.rand() < poisson.cdf(self.disease_status_days, self.DI):
                    self.status = 'J' # directly go isolation
                    self.disease_status_days = 1 # reset disease status
                    self.quarantined_days = 0 # clear quarantine days
                else:
                    self.disease_status_days += 1
                    
            # After the day of quaratine, some people may think to extend
            if self.quarantined_days > self.DQ:
                if self.infected:
                    self.status = 'E'
                else:
                    self.status = 'S'
                self.quarantined_days = 0
        
        # if isolated
        elif self.status == 'J':
            if self.disease_status_days > self.DR:
                self.status = 'R'            # being recovered 
                self.infected = False        # mark as no infected
                self.disease_status_days = 0 # reset status days
            else:
                self.disease_status_days += 1

class community:
    def __init__(self, params = default_params, structure = 'struct_free', struct_params = [1000],seed=None):
        # set default parameter
        self.params = params
        # mask rate
        self.protection_rate = params['pr']
        # paramters to generate structure
        self.struct_params = struct_params
        self.seed = seed # random seed for result
        self.lockdown_threshold = self.params['lt']
        self.lockdown_scale = self.params['ls']
        self.in_lockdown = False
        
        if structure == 'struct_free':
            #default create a community of 1000 students and 50 teachers for structure free
            self.structure = self.fully_connected(struct_params[0])
        elif structure == 'ws_network':
            #create waltz-shogatz, it require 3 parameters
            self.structure = self.watts_strogatz(self.struct_params[0],self.struct_params[1],self.struct_params[2])
        elif type(structure) == type(nx.Graph()):
            # use external defined structures
            self.structure = structure
        else:
            raise TypeError

        # initiate person nodes
        self.nodes = self.construct_nodes()

        # update protection as initial protection rate
        self.update_protection_policy()
        
    # create structure free full connected network
    def fully_connected(self, n):
        G = nx.complete_graph(n)
        # set structure label
        nx.set_node_attributes(G, -1, 'infect_source')
        nx.set_node_attributes(G, n, 'kcore')

        nx.set_node_attributes(G, 'S', 'status')
        nx.set_node_attributes(G, False, 'infected')
        nx.set_node_attributes(G, False, 'protected')
        nx.set_node_attributes(G, False, 'lockdown')
      
        return G
    
    # create watts strogatz network
    def watts_strogatz(self, n, k, p):
        G = nx.connected_watts_strogatz_graph(n, k, p, tries=100, seed=self.seed)
        #calculate k-core number
        core_num = nx.core_number(G)

        # set node attributes
        nx.set_node_attributes(G, -1, 'source')
        nx.set_node_attributes(G, core_num, 'kcore')
        nx.set_node_attributes(G, 'S', 'status')
        nx.set_node_attributes(G, False, 'infected')
        nx.set_node_attributes(G, False, 'protected')
        nx.set_node_attributes(G, False, 'lockdown')
        return G 
    
    # update protection policy about who is getting protected
    def update_protection_policy(self):
        # assign mask protection with mask rate
        for node in self.nodes:
            node.protected = (np.random.rand() <= self.protection_rate)

    # update lockdown policy
    def update_lockdown_policy(self):
        status = pd.DataFrame(self.get_all_status())
        # if the infected person reaches thresholds, impose the lockdown
        if status['status'].isin(['I','J']).mean() > self.lockdown_threshold:
            if not self.in_lockdown:
                # assign mask protection with mask rate
                for node in self.nodes:
                    node.lockdown = np.random.rand() < self.lockdown_scale
                self.in_lockdown = True
        else: # dissolve lockdown
            for node in self.nodes:
                node.lockdown = False
            self.in_lockdown = False
        return status
    
    # construct node classes
    def construct_nodes(self):
        # create node list
        node_list = [person(k, source = d['source'], kcore = d['kcore'], 
                            status =d['status'], infected = d['infected'], 
                            protected = d['protected'], params = self.params) \
                     for k, d in self.structure.nodes(data = True)]
        return node_list

    # get the clse contacs
    def get_close_contacts(self, node, status):
        return [status[i] for i in self.structure.neighbors(node)]
    
    # update for a day
    def update_status(self):
        new_nodes = self.nodes.copy()
        for i,node in enumerate(new_nodes):
            # get close contact
            close_contact = self.get_close_contacts(i, self.nodes)
            # update status
            node.update(close_contact)
        self.nodes = new_nodes
    
    # get situation of all nodes
    def get_all_status(self):
        return [{'id': node.pid, 'source': node.source, 'kcore': node.kcore, 
                 'protected': node.protected, 'status': node.status, 
                 'infected':node.infected,'spread':node.spread} \
                for node in self.nodes]
    
    # simulator of infection
    def simulator(self, days):
        count, infected, Rs = [], [], []
        com_size = float(len(self.nodes))
        for _ in range(days):
            self.update_status()
            status = self.update_lockdown_policy()
            # get status of simmulation
            count.append(status.groupby('status').size())
            infected.append(status['infected'].sum())
            Rs.append(status.loc[status['status']=='R', 'spread'].mean())
        count = pd.concat(count,axis = 1).T.fillna(0)
        #column check
        for col in ['S','E','I','J','Q','R']:
            if col not in count.columns:
                count[col] = 0
        # get numbers
        count['infected'] = infected
        count['accumulated_cases'] = count['infected'] + count['R']
        count['symptomatic_cases'] = count['I'] + count['J']
        count = count/com_size
        count['R'] = Rs
        # clean result
        count = count.fillna(0).reset_index().rename(columns ={'index':'day'})
        return count, status
        
