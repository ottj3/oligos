# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 11:36:17 2011

@author: ryan
"""

import numpy as np

class Acid:
    def __init__(self, name, p, segsize, overlapsize, minval, maxval, numpts):
        self.designtype = self.compute_design_type(numpts)
        self.locations = self.compute_locations(name, p)
        self.required_acids = self.compute_required_acids(self.locations, minval, maxval, numpts)
        acid_designs = self.compute_acid_designs(self.locations, len(p))
        segment_designs = [acidDesign_to_segmentDesign(x, len(p), segsize, overlapsize) for x in acid_designs]
        self.designlist = remove_suboptimal_designs(segment_designs)

    def compute_design_type(self, numpts):
        designtype = []
        while numpts > 1:
            numpts /= 2
            designtype.append(numpts)
        return designtype

    def compute_locations(self, name, p):
        p = list(p)
        locations = []
        while name in p:
            locations.append(p.index(name) + len(locations))
            p.remove(name)
        return locations

    def compute_required_acids(self, locations, minval, maxval, numpts):
        total_available = len(locations)
        total_required = round((maxval - minval) * total_available)
        if total_required == 0:
            total_required = 1
        delta = total_available * (maxval - minval) / (numpts - 1.0)
        points_required = [nDeltas * delta for nDeltas in self.designtype]
        required_acids = [np.floor(x) for x in points_required]
        while sum(required_acids) < total_required:
            differences = [x - y for x, y in zip(points_required, required_acids)]
            index_of_max = differences.index(max(differences))
            required_acids[index_of_max] += 1
        return [int(x) for x in required_acids]

    def compute_acid_designs(self, locations, seqlength):
        designlist = get_possible_intervals(seqlength, locations, self.required_acids[0])
        for i in range(1, len(self.required_acids)):
            designlist = disjoint_combos(designlist,
                                         get_possible_intervals(seqlength, locations, self.required_acids[i]))
        return designlist


def overlap_test(array1, array2):
    "Tests if two vectors have nonzero entries at any common index"
    return sum(array1 * array2) != 0


def disjoint_combos(list1, list2):
    """Each list consists of Design instances. Designs are tested for overlap
    and if they do not overlap are combined. All resulting combo designs are
    returned."""
    combos = []
    for x in list1:
        for y in list2:
            if overlap_test(x.score, y.score) == False:
                combos.append(x + y)
    return combos


def get_possible_intervals(seqlength, locations, number_needed):
    possible_intervals = []
    for i in range(len(locations) - number_needed + 1):
        score = np.zeros(seqlength)
        score[locations[i]:locations[i + number_needed - 1] + 1] = 1
        possible_intervals.append(Design([(locations[i], locations[i + number_needed - 1])], score))
    return possible_intervals


def acidDesign_to_segmentDesign(acid_design, seqlength, segsize, overlapsize):
    numsegs = int(np.ceil((float(seqlength) - overlapsize) / (segsize - overlapsize)))
    segment_intervals = []
    segment_score = np.ones(numsegs)
    for x, y in acid_design.intervals: #L is index of leftmost segment in design, R is index of rightmost segment
        if x <= segsize - 1:
            L = 0
        else:
            L = (x - overlapsize) / (segsize - overlapsize)
        if y >= seqlength - segsize + 1:
            R = numsegs - 1
        else:
            R = y / (segsize - overlapsize)
        segment_intervals.append((L, R))
        segment_score[L:R + 1] *= 2
    return Design(segment_intervals, segment_score)


def remove_suboptimal_designs(designlist):
    numdesigns = len(designlist)
    for i in range(numdesigns - 1, -1, -1):
        Fail = np.any([min(designlist[i].score - x.score) >= 0 and designlist[i] != x for x in designlist])
        if Fail == True:
            designlist.remove(designlist[i])
    return designlist


class Design:
    def __init__(self, intervals, score):
        self.intervals = intervals
        self.score = score

    def __add__(self, other):
        return Design(self.intervals + other.intervals, self.score + other.score)


class IndexAndScore:
    def __init__(self, index, score):
        self.index = index
        self.score = score


def compute_best_score_combo(scorelist):
    index_and_score_list = []
    for i in range(len(scorelist[0])):
        index_and_score_list.append(IndexAndScore([i, ], scorelist[0][i]))
    for j in range(1, len(scorelist)):
        index_and_score_list = compute_pairwise_score_combos(index_and_score_list, scorelist[j])
    minindex = 0
    minscorearray = index_and_score_list[0]
    minscore = sum(index_and_score_list[0].score)
    for i in range(1, len(index_and_score_list)):
        if sum(index_and_score_list[i].score) < minscore:
            minscore = sum(index_and_score_list[i].score)
            minindex = index_and_score_list[i].index
            minscorearray = index_and_score_list[i].score
    return minindex, minscore, minscorearray


def compute_pairwise_score_combos(list1, list2):
    """first argument is a list of IndexAndScore type, second list is a list of score arrays"""
    newlist = []
    for x in list1:
        for i in range(len(list2)):
            newlist.append(IndexAndScore(x.index + [i, ], x.score * list2[i]))
    return newlist


def compute_best_design(p, aoi, segsize, overlapsize, mins, maxs, numpts):
    acidlist = []
    scorelist = [[] for i in range(len(aoi))]
    for i in range(len(aoi)):
        acidlist.append(Acid(aoi[i], p, segsize, overlapsize,
                             mins[i], maxs[i], numpts[i]))
        for j in range(len(acidlist[i].designlist)):
            scorelist[i].append(acidlist[i].designlist[j].score)
    minindex, minscore, scorearray = compute_best_score_combo(scorelist)
    results = []
    if type(minindex) == int:
        results.append([aoi[0], acidlist[0].designlist[minindex].intervals, acidlist[0].required_acids])
    else:
        for i in range(len(acidlist)):
            results.append([aoi[i], acidlist[i].designlist[minindex[i]].intervals, acidlist[i].required_acids])
    return results

if __name__ == '__main__':
    print 'Initializing design.py script.'