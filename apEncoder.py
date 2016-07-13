#! /usr/bin/env python

""" Simple json 2 autoprocotol encoder. Example: cat example.json | apEncoder.py | transcriptic analyze 

Copyright (C) 2016  cweiland
Function gibsonAssembly is completely based on Scott Beckers gibsonAssembly code, see: https://github.com/scottbecker/transcriptic101/tree/master/gibson_assembly

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import os,sys,re, string, random
from pprint import pprint
import json

from autoprotocol import *
from autoprotocol.protocol import Protocol


#### custom protocol ####


def ul(microliters):
    """Unicode function name for creating microliter volumes"""
    return Unit(microliters,"microliter")                    


class CustomProtocol(Protocol):


    def __init__(self):
        #catalog inventory
        self.inv = {
            "water"       : "rs17gmh5wafm5p", # catalog; Autoclaved MilliQ H2O; ambient
            "DH5a"        : "rs16pbj944fnny", # catalog; Zymo DH5a; cold_80
            "Gibson Mix"  : "rs16pfatkggmk5", # catalog; Gibson Mix (2X); cold_20
            'NEBuilder_Master_Mix': 'rs18pc86ykcep6', # catalog; NEBuilder MasterMix
            "LB Miller"   : "rs17bafcbmyrmh", # catalog; LB Broth Miller; cold_4
            "Amp 100mgml" : "rs17msfk8ujkca", # catalog; Ampicillin 100mg/ml; cold_20
            "pUC19"       : "rs17tcqmncjfsh", # catalog; pUC19; cold_20            
            "IPTG"          : "rs18vwgfgxq597", # catalog: 100mM
            "M13_F"                       : "rs17tcpqwqcaxe", # catalog; M13 Forward (-41); cold_20 (1ul = 100pmol)
            "M13_R"                       : "rs17tcph6e2qzh", # catalog; M13 Reverse (-48); cold_20 (1ul = 100pmol)
            "SensiFAST_SYBR_No-ROX"       : "rs17knkh7526ha", # catalog; SensiFAST SYBR for qPCR      
            # kits (must be used differently)
            "lb-broth-100ug-ml-amp_6-flat" : "ki17sbb845ssx9", # catalog; ampicillin plates
            "noAB-amp_6-flat" : "ki17reefwqq3sq" # catalog; no antibiotic plates            
        }
        
        super(CustomProtocol,self).__init__()
       
    def ref(self, name, id=None, cont_type=None, storage=None, discard=None):        
        if discard:
            storage=None
        
        return super(CustomProtocol,self).ref(name, id=id, cont_type=cont_type, 
                                              storage=storage, discard=discard)
    
    
    def transfer(self, source, dest, volume, one_source=False, one_tip=False, 
                aspirate_speed=None, dispense_speed=None, 
                aspirate_source=None, dispense_target=None, 
                pre_buffer=None, disposal_vol=None, 
                transit_vol=None, blowout_buffer=None, 
                tip_type=None, new_group=False, **mix_kwargs):
        
        if type(volume) == list:
            min_volume = min(volume)
        else:
            min_volume = volume        
        
        if min_volume < ul(10) and not mix_kwargs.get('mix_after') \
           and not mix_kwargs.get('ignore_mix_after_warning'):
            raise Exception('mix_after required for <10uL of solution to ensure complete transfer. \n'
                            'Ensure you have are pipetting into something big enough and set this')
            
        if 'ignore_mix_after_warning' in mix_kwargs:
            del mix_kwargs['ignore_mix_after_warning']
            
        super().transfer(source, dest, volume, one_source=one_source, one_tip=one_tip, 
              aspirate_speed=aspirate_speed, dispense_speed=dispense_speed, 
              aspirate_source=aspirate_source, dispense_target=dispense_target, 
              pre_buffer=pre_buffer, disposal_vol=disposal_vol, 
              transit_vol=transit_vol, blowout_buffer=blowout_buffer, 
              tip_type=tip_type, new_group=new_group,**mix_kwargs)
        
    def distribute(self, source, dest, volume, allow_carryover=False,
                   mix_before=False, mix_after=False, mix_vol=None, repetitions=10,
                   flowrate="100:microliter/second", aspirate_speed=None,
                   aspirate_source=None, distribute_target=None,
                   pre_buffer=None, disposal_vol=None, transit_vol=None,
                   blowout_buffer=None, tip_type=None, new_group=False,
                   ignore_mix_after_warning=False):
        

        if volume < ul(10):
            for dest_well in dest:
                self.transfer(source,dest_well,volume,mix_before=mix_before,
                              mix_after=mix_after,mix_vol=mix_vol,
                              aspirate_speed=aspirate_speed, 
                              aspirate_source=aspirate_source, 
                              pre_buffer=pre_buffer, disposal_vol=disposal_vol, 
                              transit_vol=transit_vol, blowout_buffer=blowout_buffer,
                              tip_type=tip_type, new_group=new_group,
                              ignore_mix_after_warning=ignore_mix_after_warning)                               
        else:
            super().distribute(source, dest, volume, allow_carryover=allow_carryover,
                               mix_before=mix_before, mix_vol=mix_vol, repetitions=repetitions,
                               flowrate=flowrate, aspirate_speed=aspirate_speed,
                               aspirate_source=aspirate_source, distribute_target=distribute_target,
                               pre_buffer=pre_buffer, disposal_vol=disposal_vol, transit_vol=transit_vol,
                               blowout_buffer=blowout_buffer, tip_type=tip_type, new_group=new_group)


#### conversion tables ####

#brick to transcriptic

b2t = { 'ecoli':'rs16pbj944fnny',
        "water"       : "rs17gmh5wafm5p", # catalog; Autoclaved MilliQ H2O; ambient
        "DH5a"        : "rs16pbj944fnny", # catalog; Zymo DH5a; cold_80
        "Gibson Mix"  : "rs16pfatkggmk5", # catalog; Gibson Mix (2X); cold_20
        'NEBuilder_Master_Mix': 'rs18pc86ykcep6', # catalog; NEBuilder MasterMix
        "LB Miller"   : "rs17bafcbmyrmh", # catalog; LB Broth Miller; cold_4
        "Amp 100mgml" : "rs17msfk8ujkca", # catalog; Ampicillin 100mg/ml; cold_20
        "pUC19"       : "rs17tcqmncjfsh", # catalog; pUC19; cold_20            
        "IPTG"          : "rs18vwgfgxq597", # catalog: 100mM
        "M13_F"                       : "rs17tcpqwqcaxe", # catalog; M13 Forward (-41); cold_20 (1ul = 100pmol)
        "M13_R"                       : "rs17tcph6e2qzh", # catalog; M13 Reverse (-48); cold_20 (1ul = 100pmol)
        "SensiFAST_SYBR_No-ROX"       : "rs17knkh7526ha", # catalog; SensiFAST SYBR for qPCR      
        # kits (must be used differently)
        "lb-broth-100ug-ml-amp_6-flat" : "ki17sbb845ssx9", # catalog; ampicillin plates
        "noAB-amp_6-flat" : "ki17reefwqq3sq", # catalog; no antibiotic plates
        "EcoRI":    "rs17ta8xftpdk6",   # catalog; EcoRI-HF; cold_20; 20 units/ul
        "hindiii":    "rs18nw6kpnp44v",   # catalog; hindiii-HF; cold_20; 20 units/ul
        "CutSmart": "rs17ta93g3y85t",   # catalog; CutSmart Buffer 10x; cold_20
        #"pUC19":    "rs17tcqmncjfsh",   # catalog; pUC19; cold_20
        'reagent_plate': '', #optional: only needed after the first run              
        'NEBuilder_Master_Mix': 'rs18pc86ykcep6' # catalog; NEBuilder MasterMix
}

#### workflow switcher ####

#def workflowSwitch(n):

#### workflows ####

def dispense(protocol, species):
    plate=protocol.ref("growth_plate", id=None, cont_type="96-flat", storage="cold_4", discard=None)
    tube=protocol.ref("bacteria_tube", id=None, cont_type="micro-1.5", storage=None, discard=True)
    #species 2 tube
    s2t=protocol.provision(species, tube.well(0), "15:microliter")
    #medium 2 plate
    m2p=protocol.dispense(plate, "lb-broth-noAB", [{"column": 0, "volume": "175:microliter"}])
    return (plate,tube,s2t, m2p)


def gibsonAssembly(p, plasmid):
    
    def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
        return ''.join(random.choice(chars) for _ in range(size))
    
    experiment_name = id_generator()

    def ul(microliters):
        """Unicode function name for creating microliter volumes"""
        return Unit(microliters,"microliter")

    def expid(val,expt_name=experiment_name):
        """Generate a unique ID per experiment"""
        global experiment_name
        if not expt_name:
            assert experiment_name, "Must set experiment name"
            expt_name = experiment_name
        return "{}_{}".format(expt_name, val)

    plate_name = 're_reagents'
    number_of_experiments = 6
    
    ### step.1. provision ###
    
    dead_volume = ul(3)
    reagent_plate = p.ref(plate_name, cont_type="96-pcr", storage="cold_20")
    ecori_hindiii_well = reagent_plate.wells(["A1"])[0]

    #provisioning less than 10uL is dangerous due to pipetting mistakes (can't mix after)
    volume_per_experiment = ul(1)
    re_volume = max(ul(10),number_of_experiments*volume_per_experiment+dead_volume/2)
    p.provision(b2t["EcoRI"], ecori_hindiii_well, re_volume)
    p.provision(b2t["hindiii"], ecori_hindiii_well, re_volume)
    ecori_hindiii_well.name = 'EcoRI+hindiii'
    ecori_hindiii_well.properties = {'Unit Concentration':'10 units/ul'} 

    #there is an extra negative control for every 2 experiments so we use 3/2*experiments
    cutsmart_well = reagent_plate.wells(["B1"])[0]
    volume_per_experiment = ul(5)
    p.provision(b2t["CutSmart"], cutsmart_well, number_of_experiments*3/2*volume_per_experiment+dead_volume)
    cutsmart_well.name = 'CutSmart 10x'

    pUC19_well = reagent_plate.wells(["H12"])[0]
    volume_per_experiment = ul(1)
    pUC19_volume = max(ul(10),number_of_experiments*volume_per_experiment+dead_volume)
    p.provision(plasmid, pUC19_well, pUC19_volume)
    pUC19_well.name = 'pUC19'
    pUC19_well.properties = {'Mass Concentration':'1 ug/ul'} 

    ### step 02 linearize plasmid ###
    
    # Tubes and plates we use and then discard

    water_tube = p.ref("water_tube", cont_type="micro-1.5", discard=True).well(0)
    pcr_plate  = p.ref("pcr_plate",  cont_type="96-pcr", discard=True)

    # The result of the experiment, a pUC19 cut by EcoRI and HindIII, goes in this tube for storage
    puc19_cut_tube  = p.ref(expid("puc19_cut",experiment_name), cont_type="micro-1.5", storage="cold_20").well(0)

    # -------------------------------------------------------------
    # Provisioning and diluting.
    # Diluted EcoRI can be used more than once
    #
    p.provision(b2t["water"], water_tube, ul(500))

    # -------------------------------------------------------------
    # Restriction enzyme cutting pUC19
    # 3 experiments (1 without re)
    # 50ul total reaction volume for cutting 1ug of DNA:
    # 42uL water
    # 5ul CutSmart 10x
    # 1ul pUC19 1ml/ml or 1ug/ul (1ug of DNA)
    # 2ul EcoRI+hindiii (20 units each, >10 units per ug DNA)

    p.distribute(water_tube, pcr_plate.wells(["A1","B1","A2"]), ul(42))
    p.distribute(pUC19_well, pcr_plate.wells(["A1","B1","A2"]), ul(1),
                 mix_before=True, mix_after=True, mix_vol=ul(5))
    p.distribute(cutsmart_well, pcr_plate.wells(["A1","B1","A2"]), ul(5), 
             mix_before=True, mix_after=True, mix_vol=ul(10))
    p.distribute(water_tube, pcr_plate.wells(["A2"]), ul(2),ignore_mix_after_warning=True)
    p.distribute(ecori_hindiii_well, pcr_plate.wells(["A1","B1"]), ul(2), 
             mix_before=True, mix_after=True, mix_vol=ul(5))
    assert all(well.volume == ul(50) for well in pcr_plate.wells(["A1","B1","A2"]))

    p.mix(pcr_plate.wells(["A1","B1","A2"]), volume=ul(25), repetitions=10)

    # Incubation to induce cut, then heat inactivation of EcoRI+hindiii
    p.seal(pcr_plate)
    #NEB enzymes are time saver so they only need 15minutes
    p.incubate(pcr_plate, "warm_37", "15:minute", shaking=False)
    #inactivate RE's (see this site for thermal temps)
    #https://www.neb.com/tools-and-resources/usage-guidelines/heat-inactivation
    p.thermocycle(pcr_plate, [{"cycles":  1, "steps": [{"temperature": "80:celsius", "duration": "20:minute"}]}], volume=ul(50))

    # --------------------------------------------------------------
    # Gel electrophoresis, to ensure the cutting worked
    #
    p.unseal(pcr_plate)
    p.mix(pcr_plate.wells(["A1","B1","A2"]), volume=ul(25), repetitions=5)
    p.transfer(water_tube, pcr_plate.wells(["D1","E1","D2"]), ul(15))
    p.transfer(pcr_plate.wells(["A1","B1","A2"]), pcr_plate.wells(["D1","E1","D2"]), ul(8), mix_after=True, mix_vol=ul(10))
    assert all(well.volume == ul(20) + dead_volume for well in pcr_plate.wells(["D1","E1","D2"]))
    
    p.gel_separate(pcr_plate.wells(["D1","E1","D2"]), ul(20), "agarose(10,2%)", "ladder2", "15:minute", expid("gel",experiment_name))
    remaining_volumes = [well.volume - dead_volume for well in pcr_plate.wells(["A1","B1"])]
    # Then consolidate all cut plasmid to one tube (puc19_cut_tube).
    p.consolidate(pcr_plate.wells(["A1","B1"]), puc19_cut_tube, remaining_volumes, allow_carryover=True)

    #verify that we haven't gone below what can actually be pulled from tubes
    #assert_valid_volume([water_tube, ecori_hindiii_well, cutsmart_well, pUC19_well])

    ### step 03 full Gibson assembly and transformation protocol for sfGFP and pUC19 ###

    sfgfp_pcroe_amp_tube = p.ref("pfragment", cont_type="micro-1.5", storage="cold_20").well(0)    
    clone_plate = p.ref(expid("clone"), cont_type="96-pcr", storage="cold_20")
    water_tube = p.ref("water", cont_type="micro-1.5", discard=True).well(0)

    ### 03.01 do transformate ###     
    
    #
    # Combine all the Gibson reagents in one tube and thermocycle
    #

    p.provision(p.inv["NEBuilder_Master_Mix"], clone_plate.well(0), ul(10))
    p.transfer(water_tube, clone_plate.well(0), ul(3.5), mix_after=True, mix_vol=ul(6))
    p.transfer(puc19_cut_tube, clone_plate.well(0), ul(2.5), mix_after=True, mix_vol=ul(6))
    p.transfer(sfgfp_pcroe_amp_tube, clone_plate.well(0), ul(4), mix_after=True, mix_vol=ul(10))

    p.seal(clone_plate)
    p.thermocycle(clone_plate, [{"cycles":  1, "steps": [{"temperature": "50:celsius", "duration": "15:minute"}]}], \
                                                     volume=ul(20))
    #
    # Dilute assembled plasmid 4X according to the NEB Gibson assembly protocol (20ul->80ul)
    #
    p.unseal(clone_plate)
    p.transfer(water_tube, clone_plate.well(0), ul(60), mix_after=True, mix_vol=ul(40), repetitions=5)

                       
    
data = json.load(sys.stdin)
data.sort(key=lambda x:x['id'])

(species, method) = (None,None)
p = CustomProtocol()

for item in data:
    if item['elementtype'] == 'organism' or item['elementtype'] == 'plasmid':
        brick = item['name']
        # lookup catalog for tsc reference
        if brick in b2t.keys():
            brick = b2t[brick]
    elif item['elementtype'] == 'bioprotocol':
        method = item['name']
    else:
        sys.exit(0)

    #trigger if brick and method    
    if brick and method:    
        exec("%s(%s,\"%s\")" % (method,"p",brick))
        #print("%s(%s,\"%s\")" % (method,"p",brick)) 
# Dump autoprotocol2JSON.
p2ap = json.dumps(p.as_dict(), indent=2)
print(p2ap)

