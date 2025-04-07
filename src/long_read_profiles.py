############################################################################
# Copyright (c) 2022-2024 University of Helsinki
# Copyright (c) 2020-2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
from functools import partial
from collections import defaultdict
from itertools import combinations

from .common import (
    contains,
    contains_approx,
    equal_ranges,
    junctions_from_blocks,
    left_of,
    overlaps,
    overlaps_at_least,
    overlaps_at_least_when_overlap,
    interval_bin_search,
    interval_bin_search_rev
)

logger = logging.getLogger('IsoQuant')


class MappedReadProfile:
    def __init__(self, gene_profile, read_profile, read_features, gene_profile_range=None):
        self.gene_profile = gene_profile
        self.read_profile = read_profile
        self.read_features = read_features
        # profile range points to indices in the mapped region
        self.gene_profile_range = (0, len(self.gene_profile)) if gene_profile_range is None else gene_profile_range


class CombinedReadProfiles:
    def __init__(self, read_intron_profile, read_exon_profile, read_split_exon_profile,
                 polya_info=None, cage_hits=-1, unique_imputation=True, alignment=None):
        self.read_intron_profile = read_intron_profile
        self.read_exon_profile = read_exon_profile
        self.read_split_exon_profile = read_split_exon_profile
        self.alignment = alignment
        self.polya_info = polya_info
        self.cage_hits = cage_hits
        self.unique_imputation = unique_imputation
## This class is used to impute exon intervals when coverage information is missing in BaseCode data.
class ExonImputation:
    def __init__(self, known_features, gene_region,
                 comparator = contains_approx,
                 overlap = overlaps,
                 delta=0):
        self.known_features = known_features
        self.gene_region = gene_region
        self.comparator = comparator
        self.overlap = overlap
        self.delta = delta
    def overlaps(self, range1, range2):
        return not (range1[1] < range2[0] or range1[0] > range2[1])
    def invert_introns(self, matches):
        new_end = matches[0][0] - 1
        new_start = matches[-1][1] + 1
        tmp_new_blocks = []
        for i in range(len(matches)-1):
            tmp_new_blocks.append((matches[i][1]+1, matches[i+1][0]-1))
        return new_end, tmp_new_blocks, new_start
    def prune_exons_outside_gene_model(self,sorted_blocks, sorted_deleted_blocks):
        prune_left = 0
        prune_right = 0
        n_blocks = len(sorted_blocks)
        n_blocks_deleted = len(sorted_deleted_blocks)
        for i in range(n_blocks-1):
            if i >= n_blocks_deleted:
                break
            current_block_exon = sorted_blocks[i]
            if not overlaps(current_block_exon, self.gene_region):
                if current_block_exon[1] == sorted_deleted_blocks[i][0] - 1:
                    prune_left += 1
                else:
                    break
            else:
                break
        for i in range(1,n_blocks):
            current_block_exon = sorted_blocks[-i]
            if i >= n_blocks_deleted:
                break
            if not self.overlaps(current_block_exon, self.gene_region):
                if current_block_exon[0] == sorted_deleted_blocks[-i][1] + 1:
                    prune_right += 1
                else:
                    break
            else:
                break
        return sorted_blocks[prune_left:(n_blocks-prune_right)], sorted_deleted_blocks[prune_left:(n_blocks_deleted-prune_right)]
    def impute_exon_structure(self, sorted_blocks, sorted_deleted_blocks):
        ### TODO: Validate and test function
        if overlaps((sorted_blocks[0][0], sorted_blocks[-1][1]), self.gene_region):
            sorted_blocks, sorted_deleted_blocks = self.prune_exons_outside_gene_model(sorted_blocks, sorted_deleted_blocks)
        #print(sorted_blocks)
        if not sorted_deleted_blocks:
            #print('no deletions')
            return sorted_blocks, True
        new_sorted_blocks = []
        unique_imputation = True
        exon_pos = 0
        deletion_pos = 0
        new_start = None
        start = -1
        while exon_pos < len(sorted_blocks) and deletion_pos < len(sorted_deleted_blocks):
            # print('Exon pos: ', exon_pos)
            # print('Deletion pos: ', deletion_pos)
            # print(start)
            current_block = sorted_blocks[exon_pos]
            if new_start == -1 or start == -1:
                # print('Changing start to', current_block[0])
                start = current_block[0]
            deleted_block = sorted_deleted_blocks[deletion_pos]
            if current_block[0] < deleted_block[0] or (deletion_pos + 1) == len(sorted_deleted_blocks):
                # The current block is before the next deleted block or there are no more deleted blocks
                if (sorted_blocks[exon_pos][1] + 1) == sorted_deleted_blocks[deletion_pos][0]:
                    # The exon block is followed by a deletion
                    deleted_block = sorted_deleted_blocks[deletion_pos]
                    new_end, tmp_new_blocks, new_start, unique_imputation = self.find_features_in_block(deleted_block, current_block)
                    # print(new_end, tmp_new_blocks, new_start, unique_imputation)
                    if new_start is not None:
                        # print('Adding block from deletion: ', (start, new_end))
                        new_sorted_blocks.append((start, new_end))
                        new_sorted_blocks.extend(tmp_new_blocks)
                        start = new_start
                else:
                    # The exon block is followed by a refskip or is the last block
                    end = current_block[1]
                    # print('Adding block from refskip or last block: ', (start, end))
                    new_sorted_blocks.append((start, end))
                    new_start = -1
                exon_pos += 1
            else:
                # The current block is after the next deleted block
                deletion_pos += 1
        return new_sorted_blocks, unique_imputation
    def find_features_in_block(self, deleted_block, current_block):
        # print('Deleted block: ', deleted_block)
        unique_imputation = True
        matches = 0
        match_list = []
        for gene_pos in range(len(self.known_features)):
            isoform_feature = self.known_features[gene_pos]
            # print('Intron: ', isoform_feature)
            if self.comparator(deleted_block, isoform_feature):
                # print('Match: ', isoform_feature)
                matches += 1
                match_list.append(isoform_feature)
        if matches == 0:
            for gene_pos in range(len(self.known_features)):
                isoform_feature = self.known_features[gene_pos]
            # print('Intron: ', isoform_feature)
                if self.overlap(deleted_block, isoform_feature):
                    unique_imputation = False
            # No intron in deleted block
            new_end = None 
            tmp_new_blocks = []
            new_start = None
        elif matches == 1:
            # One intron found in the deleted block
            new_end = match_list[0][0] - 1
            tmp_new_blocks = []
            new_start = match_list[0][1] + 1
        else:
            # Check if there is any overlap between contained introns
            if any(overlaps(intron1, intron2) for intron1, intron2 in combinations(match_list, 2)):
                new_end = min([t[0] - 1 for t in match_list])
                tmp_new_blocks = []
                new_start = max([t[1] + 1 for t in match_list])
                unique_imputation = False
            else:
                new_end, tmp_new_blocks, new_start = self.invert_introns(match_list)
        return new_end, tmp_new_blocks, new_start, unique_imputation

# The following 2 classes are very similar, but lets keep them separately for now
# accepts sorted gapless alignment blocks
class OverlappingFeaturesProfileConstructor:
    # ignore_terminal -- bool flag, indicates whether to ignore leading and trailing -1s in the profile
    def __init__(self, known_features, gene_region,
                 comparator = partial(equal_ranges, delta=0),
                 absence_condition = contains,
                 delta=0):
        self.known_features = known_features
        self.gene_region = gene_region
        self.comparator = comparator
        self.absence_condition = absence_condition
        self.delta = delta
    
    

    def construct_intron_profile(self, sorted_blocks, polya_position=-1, polyt_position=-1):
        if not sorted_blocks:
            return  MappedReadProfile([], [], defaultdict(list), (0, 0))
        mapped_region = (sorted_blocks[0][0], sorted_blocks[-1][1])
        read_introns = junctions_from_blocks(sorted_blocks)
        return self.construct_profile_for_features(read_introns, mapped_region, polya_position, polyt_position)

    def construct_exon_profile(self, sorted_blocks, polya_position=-1, polyt_position=-1):
        if not sorted_blocks:
            return  MappedReadProfile([], [], defaultdict(list), (0, 0))
        mapped_region = (sorted_blocks[0][1] + self.delta, sorted_blocks[-1][0] - self.delta)
        return self.construct_profile_for_features(sorted_blocks, mapped_region, polya_position, polyt_position)

    def match_delta(self, feature1, feature2):
        return abs(feature1[0] - feature2[0]) + abs(feature1[1] - feature2[1])

    def match_genomic_features(self, read_features):
        matched_features = defaultdict(list)

        # TODO: starting value can be detected using binary search for long profiles
        gene_pos = 0
        read_pos = 0
        while gene_pos < len(self.known_features) and read_pos < len(read_features):
            if self.comparator(read_features[read_pos], self.known_features[gene_pos]):
                matched_features[read_pos].append(gene_pos)
                gene_pos += 1
            elif overlaps(read_features[read_pos], self.known_features[gene_pos]):
                gene_pos += 1
            elif left_of(read_features[read_pos], self.known_features[gene_pos]):
                read_pos += 1
            else:
                gene_pos += 1
        # eliminating non unique features
        for read_pos in matched_features.keys():
            if len(matched_features[read_pos]) > 1:
                deltas = [self.match_delta(read_features[read_pos], self.known_features[gene_pos])
                          for gene_pos in matched_features[read_pos]]
                best_match = min(deltas)
                filtered_features = []
                for i, d in enumerate(deltas):
                    if d == best_match:
                        filtered_features.append(matched_features[read_pos][i])
                matched_features[read_pos] = filtered_features

        corrected_features = []
        for i, intron in enumerate(read_features):
            if i in matched_features:
                # known feature
                # corresponding known feature, always take first for now
                corrected_features.append(self.known_features[matched_features[i][0]])
            else:
                corrected_features.append(intron)
        return corrected_features

    def construct_profile_for_features(self, read_features, mapped_region=(0, 0), polya_position=-1, polyt_position=-1):
        read_profile = [0] * (len(read_features))
        intron_profile = [0] * (len(self.known_features))
        matched_features = defaultdict(list)

        for i in range(len(intron_profile)):
            if self.absence_condition(mapped_region, self.known_features[i]):
                intron_profile[i] = -1
        for i in range(len(read_profile)):
            if self.absence_condition(self.gene_region, read_features[i]):
                read_profile[i] = -1

        gene_pos = 0
        read_pos = 0
        while gene_pos < len(self.known_features) and read_pos < len(read_features):
            read_feature = read_features[read_pos]
            isoform_feature = self.known_features[gene_pos]
            if read_feature[1] < isoform_feature[0]:
                if read_profile[read_pos] == 0 and gene_pos > 0:
                    read_profile[read_pos] = -1
                read_pos += 1
            elif isoform_feature[1] < read_feature[0]:
                if read_pos > 0:
                    intron_profile[gene_pos] = -1
                gene_pos += 1
            elif self.comparator(read_feature, isoform_feature):
                intron_profile[gene_pos] = 1
                read_profile[read_pos] = 1
                matched_features[read_pos].append(gene_pos)
                gene_pos += 1
            elif overlaps(read_feature, isoform_feature):
                if self.absence_condition(mapped_region, isoform_feature):
                    intron_profile[gene_pos] = -1
                gene_pos += 1

        # eliminating non unique features
        for read_pos in matched_features.keys():
            if len(matched_features[read_pos]) > 1:
                deltas = [self.match_delta(read_features[read_pos], self.known_features[gene_pos])
                          for gene_pos in matched_features[read_pos]]
                best_match = min(deltas)
                for i in range(len(matched_features[read_pos])):
                    if deltas[i] > best_match:
                        intron_profile[matched_features[read_pos][i]] = -1

        corrected_features = []
        for i, v in enumerate(read_profile):
            if v == 1:
                # known feature
                # corresponding known feature, always take first for now
                corrected_features.append(self.known_features[matched_features[i][0]])
            else:
                corrected_features.append(read_features[i])

        # making everying beyond polyA tail as outside feature
        if polya_position != -1:
            for i in range(len(self.known_features)):
                # feature is surely beyond polyA tail
                if self.known_features[i][0] > polya_position + self.delta:
                    intron_profile[i] = -2
        if polyt_position != -1:
            for i in range(len(self.known_features)):
                # feature is surely before polyT tail
                if self.known_features[i][1] < polyt_position - self.delta:
                    intron_profile[i] = -2

        start_pos = 0
        while start_pos < len(intron_profile) and intron_profile[start_pos] == 0:
            start_pos += 1
        end_pos = len(intron_profile) - 1
        while end_pos >= 0 and intron_profile[end_pos] == 0:
            end_pos -= 1

        return MappedReadProfile(intron_profile, read_profile, read_features, (start_pos, end_pos + 1))


# accepts sorted gapless alignment blocks
class NonOverlappingFeaturesProfileConstructor:
    def __init__(self, known_exons, comparator=overlaps, delta=0):
        self.known_exons = known_exons
        self.comparator = comparator
        self.delta = delta

    def construct_profile(self, sorted_blocks, polya_position=-1, polyt_position=-1):
        if not sorted_blocks:
            return  MappedReadProfile([], [], defaultdict(list), (0, 0))
        exon_profile = [0] * (len(self.known_exons))
        read_profile = [0] * (len(sorted_blocks))
        read_exons = sorted_blocks
        gene_pos = 0
        read_pos = 0

        while gene_pos < len(self.known_exons) and read_pos < len(read_exons):
            read_exon = read_exons[read_pos]
            gene_exon = self.known_exons[gene_pos]
            if read_exon[1] < gene_exon[0]:
                if gene_pos > 0 and read_profile[read_pos] == 0:
                    read_profile[read_pos] = -1
                read_pos += 1
            elif gene_exon[1] < read_exon[0]:
                if read_pos > 0 and exon_profile[gene_pos] == 0:
                    exon_profile[gene_pos] = -1
                gene_pos += 1
            elif self.comparator(read_exon, gene_exon):
                exon_profile[gene_pos] = 1
                read_profile[read_pos] = 1
                if read_exon[1] < gene_exon[1]:
                    read_pos += 1
                else:
                    gene_pos += 1
            else:
                if read_exon[1] < gene_exon[1]:
                    read_pos += 1
                else:
                    gene_pos += 1

        # making everything beyond polyA tail as outside feature
        if polya_position != -1:
            polya_index = interval_bin_search(self.known_exons, polya_position + self.delta)
            if polya_index != -1:
                for i in range(polya_index + 1, len(self.known_exons)):
                    exon_profile[i] = -2
        if polyt_position != -1:
            polyt_index = interval_bin_search_rev(self.known_exons, polyt_position - self.delta)
            if polyt_index != -1:
                for i in range(0, polyt_index):
                    exon_profile[i] = -2

        start_pos = 0
        while start_pos < len(exon_profile) and exon_profile[start_pos] == 0:
            start_pos += 1
        end_pos = len(exon_profile) - 1
        while end_pos >= 0 and exon_profile[end_pos] == 0:
            end_pos -= 1

        return MappedReadProfile(exon_profile, read_profile, read_exons, (start_pos, end_pos + 1))


class CombinedProfileConstructor:
    def __init__(self, gene_info, params):
        self.gene_info = gene_info
        self.params = params
        gene_region = (gene_info.start, gene_info.end)

        self.exon_imputation = ExonImputation(self.gene_info.intron_profiles.features, gene_region, comparator=contains, delta=self.params.delta)

        self.intron_profile_constructor = \
            OverlappingFeaturesProfileConstructor(self.gene_info.intron_profiles.features, gene_region,
                                                  comparator=partial(equal_ranges, delta=self.params.delta),
                                                  absence_condition=partial(overlaps_at_least, delta=self.params.minimal_intron_absence_overlap),
                                                  delta=self.params.delta)
        self.exon_profile_constructor = \
            OverlappingFeaturesProfileConstructor(self.gene_info.exon_profiles.features, gene_region,
                                                  comparator=partial(equal_ranges, delta=self.params.delta),
                                                  delta=self.params.delta)
        self.split_exon_profile_constructor = \
            NonOverlappingFeaturesProfileConstructor(self.gene_info.split_exon_profiles.features,
                                                     comparator=partial(overlaps_at_least_when_overlap,
                                                                        delta=self.params.minimal_exon_overlap),
                                                     delta=self.params.delta)

    def construct_profiles(self, sorted_blocks, sorted_deleted_blocks, polya_info, cage_hits):
        #print('Before: ', sorted_blocks, sorted_deleted_blocks)
        #print(self.exon_imputation.gene_region)
        imputed_sorted_blocks, unique_imputation = self.exon_imputation.impute_exon_structure(sorted_blocks, sorted_deleted_blocks)
        if not imputed_sorted_blocks:
            print("Before: ", sorted_blocks, sorted_deleted_blocks, "After: ", imputed_sorted_blocks, "Gene Region: {}-{}".format(self.exon_imputation.gene_region[0], self.exon_imputation.gene_region[1]) )
        #print('After: ', imputed_sorted_blocks)
        intron_profile = self.intron_profile_constructor.construct_intron_profile(imputed_sorted_blocks,
                                                                                  polya_info.external_polya_pos,
                                                                                  polya_info.external_polyt_pos)
        exon_profile = None
        if self.params.count_exons:
            exon_profile = self.exon_profile_constructor.construct_exon_profile(imputed_sorted_blocks,
                                                                                polya_info.external_polya_pos,
                                                                                polya_info.external_polyt_pos)
        split_exon_profile = self.split_exon_profile_constructor.construct_profile(imputed_sorted_blocks,
                                                                                   polya_info.external_polya_pos,
                                                                                   polya_info.external_polyt_pos)
        return imputed_sorted_blocks, unique_imputation, CombinedReadProfiles(intron_profile, exon_profile, split_exon_profile,
                                    polya_info=polya_info, cage_hits=cage_hits, unique_imputation=unique_imputation)