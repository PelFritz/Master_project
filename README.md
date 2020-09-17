# Impact of genetic variation onactivity of plant gene promoters innatural accessions and cropcultivars
## Deep Learning meets plant genomics

This repository contains the scripts for my master project known in german as the "Masterarbeit".
*scenario_1.py* and *scenario_4.py* correspond to the screnarios described by the project. The other two scenarios left out can
easily be reproduce by tweaking the code a bit. For example the piece of code directly below will reproduce *scenario_2*.

`if re.search(hexamer, sequence) and re.search(hexamer2, sequence):
    prom_seq.append(sequence)
    label.append(1)
    mut_sequence = re.sub(hexamer, mut_hexamer, sequence, count=1)
    mut_sequence_2 = re.sub(hexamer2, mut_hexamer_2, mut_sequence, count=1)
    prom_seq.append(mut_sequence_2)
    label.append(0)`
