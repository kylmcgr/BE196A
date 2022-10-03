# Complete the rest of the staples in scadnano
# ALSO ENDED UP CHANGING THE EYES LATER ON


import scadnano as sc
import math

def main() -> None:
    design = create_design()
    design.write_scadnano_file()


def create_design():
    helices = [sc.Helix(max_offset=224) for _ in range(30)]
    design = sc.Design(helices=helices, grid=sc.square)

    add_scaffold_precursors(design)
    add_scaffold_crossovers(design)

    add_staple_precursors(design)
    add_staple_crossovers(design)

    return design

def add_scaffold_precursors(design: sc.Design) -> None:
    for helix in range(0, 5, 2): # top face
        design.strand(helix, 0).move(224).as_scaffold()
        design.strand(helix + 1, 224).move(-224).as_scaffold()
    design.strand(6, 0).move(224).as_scaffold()

    for helix in range(7, 12, 2): # eyes
        design.strand(helix, 59).move(-59).as_scaffold()
        design.strand(helix + 1, 0).move(59).as_scaffold()
        design.strand(helix, 144).move(-64).as_scaffold()
        design.strand(helix + 1, 80).move(64).as_scaffold()
        design.strand(helix, 224).move(-59).as_scaffold()
        design.strand(helix + 1, 165).move(59).as_scaffold()

    for helix in range(13, 16, 2): # inside face
        design.strand(helix, 224).move(-224).as_scaffold()
        design.strand(helix + 1, 0).move(224).as_scaffold()

    for helix in range(17, 24, 2): # mouth outside parts
        design.strand(helix, (27+math.floor((helix-17)/2*(10+2/3)))).move(-(27+math.floor((helix-17)/2*(10+2/3)))).as_scaffold()
        design.strand(helix + 1, 0).move((27+math.floor((helix-17)/2*(10+2/3)))).as_scaffold()
        design.strand(helix, 224).move(-(224-(197-math.floor((helix-17)/2*(10+2/3))))).as_scaffold()
        design.strand(helix + 1, (197-math.floor((helix-17)/2*(10+2/3)))).move(224-(197-math.floor((helix-17)/2*(10+2/3)))).as_scaffold()

     # mouth inside parts
    design.strand(17, 155).move(-86).as_scaffold()
    design.strand(18, 69).move(86).as_scaffold()

    for helix in range(19, 22, 2): # mouth inside parts
        left = 59+math.floor((helix-19)/2*(10+2/3))
        right = 165-math.floor((helix-19)/2*(10+2/3))
        diff = right-left
        design.strand(helix, right).move(-diff).as_scaffold()
        design.strand(helix + 1, left).move(diff).as_scaffold()

    for helix in range(25, 28, 2): # bottom face
        design.strand(helix, 224).move(-224).as_scaffold()
        design.strand(helix + 1, 0).move(224).as_scaffold()

    design.strand(29, 224).move(-112).as_scaffold() # bottom part of scaffold has a "nick"
    design.strand(29, 112).move(-112).as_scaffold()

def add_scaffold_crossovers(design: sc.Design) -> None:
    for helix in range(1, 4, 2):
        design.add_full_crossover(helix=helix, helix2=helix + 1, offset=112, forward=False) # top middle

    design.add_full_crossover(helix=5, helix2=6, offset=165, forward=False) # top right ish

    for helix in range(7, 23, 2):
        design.add_full_crossover(helix=helix, helix2=helix - 1, offset=96, forward=False) # middle of face

    for helix in range(7, 12, 2):
        design.add_half_crossover(helix=helix, helix2=helix + 1, offset=58, forward=False) # left eye left part
        design.add_half_crossover(helix=helix, helix2=helix + 1, offset=165, forward=False) # right eye right part
        design.add_half_crossover(helix=helix, helix2=helix + 1, offset=80, forward=False) # left eye right part
        design.add_half_crossover(helix=helix, helix2=helix + 1, offset=143, forward=False) # right eye left part

    for helix in range(13, 16, 2):
        design.add_full_crossover(helix=helix, helix2=helix + 1, offset=58, forward=False) # left eye left part
        design.add_full_crossover(helix=helix, helix2=helix + 1, offset=165, forward=False) # right eye right part

    for helix in range(17, 24, 2):
        design.add_half_crossover(helix=helix, helix2=helix + 1, offset=26+math.floor((helix-17)/2*(10+2/3)), forward=False) # mouth left part
        design.add_half_crossover(helix=helix, helix2=helix + 1, offset=197-math.floor((helix-17)/2*(10+2/3)), forward=False) # mouth right part

    design.add_half_crossover(helix=17, helix2=18, offset=69, forward=False) # mouth left inside part
    design.add_half_crossover(helix=17, helix2=18, offset=154, forward=False) # mouth right inside part

    for helix in range(19, 22, 2):
        design.add_half_crossover(helix=helix, helix2=helix + 1, offset=59+math.floor((helix-19)/2*(10+2/3)), forward=False) # mouth left inside part
        design.add_half_crossover(helix=helix, helix2=helix + 1, offset=164-math.floor((helix-19)/2*(10+2/3)), forward=False) # mouth right inside part

    for helix in range(25, 28, 2):
        design.add_full_crossover(helix=helix, helix2=helix + 1, offset=112, forward=False) # bottom middle

    for helix in range(0, 29, 2):  # scaffold edges crossovers
        design.add_half_crossover(helix=helix, helix2=helix + 1, offset=0, forward=True)
        design.add_half_crossover(helix=helix, helix2=helix + 1, offset=223, forward=True)

def add_staple_precursors(design: sc.Design) -> None:
    staples = [sc.Strand([sc.Domain(helix=helix, forward=helix % 2 == 1, start=0, end=224)]) for helix in range(7)] # top face

    staples += [sc.Strand([sc.Domain(helix=helix, forward=helix % 2 == 1, start=0, end=59)]) for helix in range(7, 13)] # eyes
    staples += [sc.Strand([sc.Domain(helix=helix, forward=helix % 2 == 1, start=80, end=144)]) for helix in range(7, 13)] # eyes
    staples += [sc.Strand([sc.Domain(helix=helix, forward=helix % 2 == 1, start=165, end=224)]) for helix in range(7, 13)] # eyes

    staples += [sc.Strand([sc.Domain(helix=helix, forward=helix % 2 == 1, start=0, end=224)]) for helix in range(13, 17)] # inside face

    staples += [sc.Strand([sc.Domain(helix=helix, forward=helix % 2 == 1, start=0, end=27+math.floor((helix-17-((helix+1)%2))/2*(10+2/3)))]) for helix in range(17, 25)] # mouth
    staples += [sc.Strand([sc.Domain(helix=helix, forward=helix % 2 == 1, start=197-math.floor((helix-17-((helix+1)%2))/2*(10+2/3)), end=224)]) for helix in range(17, 25)] # mouth
    staples += [sc.Strand([sc.Domain(helix=helix, forward=helix % 2 == 1, start=59+math.floor((helix-19-((helix+1)%2))/2*(10+2/3)), end=165-math.floor((helix-19-((helix+1)%2))/2*(10+2/3)))]) for helix in range(19, 23)] # mouth
    staples += [sc.Strand([sc.Domain(helix=helix, forward=helix % 2 == 1, start=69, end=155)]) for helix in range(17, 19
    )]

    staples += [sc.Strand([sc.Domain(helix=helix, forward=helix % 2 == 1, start=0, end=224)]) for helix in range(25, 30)] # bottom face

    for staple in staples:
        design.add_strand(staple)

def add_staple_crossovers(design: sc.Design) -> None:
    for helix in range(0, 29, 2):
        design.add_full_crossover(helix=helix, helix2=helix + 1, offset=16, forward=False)
        design.add_full_crossover(helix=helix, helix2=helix + 1, offset=208, forward=False)

    for helix in range(0, 6, 2):
        for i in range(2):
            design.add_full_crossover(helix=helix, helix2=helix + 1, offset=48+32*i, forward=False)
            design.add_full_crossover(helix=helix, helix2=helix + 1, offset=176-32*i, forward=False)
        for i in range(6):
            design.add_full_crossover(helix=helix + 1, helix2=helix + 2, offset=32+32*i, forward=True)

    for helix in range(25, 29, 2):
        for i in range(2):
            design.add_full_crossover(helix=helix + 1, helix2=helix + 2, offset=48+32*i, forward=False)
            design.add_full_crossover(helix=helix + 1, helix2=helix + 2, offset=176-32*i, forward=False)
        for i in range(6):
            design.add_full_crossover(helix=helix, helix2=helix + 1, offset=32+32*i, forward=True)


if __name__ == '__main__':
    main()
