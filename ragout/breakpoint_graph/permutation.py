"""
This module provides PermutationContainer class
which stores permutations and provides some filtering
procedures
"""

from collections import defaultdict
import logging
import os
import math

from ragout.shared.debug import DebugConfig

logger = logging.getLogger()
debugger = DebugConfig.get_instance()


class PermException(Exception):
    pass

class Block:
    def __init__(self, block_id, sign, start=None, end=None):
        self.block_id = block_id
        self.sign = sign
        self.start = start
        self.end = end

    def length(self):
        if not self.start or self.end:
            return None

        assert self.end >= self.start
        return self.end - self.start

    def signed_id(self):
        return self.block_id * self.sign


class Permutation:
    def __init__(self, genome_name, chr_name, chr_id, chr_len, blocks):
        self.genome_name = genome_name
        self.chr_name = chr_name
        self.chr_id = chr_id
        self.chr_len = chr_len
        self.blocks = blocks

    def iter_pairs(self):
        for pb, nb in zip(self.blocks[:-1], self.blocks[1:]):
            yield pb, nb


class PermutationContainer:
    def __init__(self, block_coords_file, recipe):
        """
        Parses permutation files referenced from recipe and filters duplications
        """
        self.ref_perms = []
        self.target_perms = []

        logging.info("Reading permutation file")
        permutations = _parse_blocks_coords(block_coords_file)
        if not permutations:
            raise PermException("Error reading permutations")

        for p in permutations:
            if p.genome_name not in recipe["genomes"]:
                continue
            if p.genome_name == recipe["target"]:
                self.target_perms.append(p)
            else:
                self.ref_perms.append(p)

        logger.debug("Read {0} reference sequences"
                     .format(len(self.ref_perms)))
        if not len(self.ref_perms):
            raise PermException("No synteny blocks found in "
                                "reference sequences")
        logger.debug("Read {0} target sequences"
                     .format(len(self.target_perms)))
        if not len(self.target_perms):
            raise PermException("No synteny blocks found in "
                                "target sequences")

        self.target_blocks = set()
        for perm in self.target_perms:
            self.target_blocks |= set(map(lambda b: b.block_id, perm.blocks))

        #filter dupilcated blocks
        self.duplications = _find_duplications(self.ref_perms,
                                               self.target_perms)
        to_hold = self.target_blocks - self.duplications
        self.ref_perms_filtered = [_filter_perm(p, to_hold)
                                      for p in self.ref_perms]
        self.target_perms_filtered = [_filter_perm(p, to_hold)
                                         for p in self.target_perms]
        self.target_perms_filtered = list(filter(lambda p: p.blocks,
                                                 self.target_perms_filtered))

        if debugger.debugging:
            file = os.path.join(debugger.debug_dir, "used_contigs.txt")
            _write_permutations(self.target_perms_filtered, open(file, "w"))


def _find_duplications(ref_perms, target_perms):
    """
    Find duplicated blocks
    """
    index = defaultdict(set)
    duplications = set()
    for perm in ref_perms + target_perms:
        for block in perm.blocks:
            if perm.genome_name in index[block.block_id]:
                duplications.add(block.block_id)
            else:
                index[block.block_id].add(perm.genome_name)

    return duplications


def _filter_perm(perm, to_hold):
    """
    Filters duplications
    """
    new_perm = Permutation(perm.genome_name, perm.chr_name, perm.chr_id,
                           perm.chr_len, [])
    for block in perm.blocks:
        if block.block_id in to_hold:
            new_perm.blocks.append(block)
    return new_perm


def _parse_permutations(filename):
    """
    Parses a file with signed permutations
    """
    permutations = []
    chr_count = 0
    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                tokens = line[1:].split(".", 1)
                if len(tokens) != 2:
                    PermException("permutation ids in " + filename + " do not "
                                 "follow naming convention: 'genome.chromosome'")
                    return None

                genome_name, chr_name = tokens
            else:
                permutations.append(Permutation(genome_name, chr_name,
                                    chr_count, None, []))
                blocks_ids = map(int, line.split(" ")[:-1])
                for b in blocks_ids:
                    permutations[-1].blocks.append(Block(abs(b),
                                                         math.copysign(1, b)))
                chr_count += 1
    return permutations



def _parse_blocks_coords(filename):
    """
    Parses a file with blocks coords
    """
    perm_by_id = {}
    with open(filename, "r") as f:
        header = True
        for line in f:
            line = line.strip()
            if not line:
                continue

            if header:
                if line.startswith("Seq_id"):
                    continue

                if line.startswith("-"):
                    header = False
                    continue

                chr_id, chr_size, seq_name = line.split("\t")
                tokens = seq_name.split(".", 1)
                if len(tokens) != 2:
                    PermException("permutation ids in " + filename + " do not "
                                  "follow naming convention: 'genome.chromosome'")
                    return None
                genome_name, chr_name = tokens
                perm_by_id[chr_id] = Permutation(genome_name, chr_name,
                                                 chr_id, int(chr_size), [])

            else:
                if line.startswith("Seq_id") or line.startswith("-"):
                    continue

                if line.startswith("Block"):
                    block_id = int(line.split(" ")[1][1:])
                    continue

                seq_id, sign, start, end, length = line.split("\t")
                if sign == "-":
                    start, end = end, start
                if int(end) < int(start):
                    raise PermException("Error in permutations file format")

                sign_num = 1 if sign == "+" else -1
                perm_by_id[seq_id].blocks.append(Block(block_id, sign_num,
                                                      int(start), int(end)))

    for perm in perm_by_id.values():
        perm.blocks.sort(key=lambda b: b.start)

    out_perms = list(filter(lambda b: len(b.blocks), perm_by_id.values()))
    return out_perms


def _write_permutations(permutations, out_stream):
    """
    Outputs permutations
    """
    for perm in permutations:
        out_stream.write(">" + perm.chr_name + "\n")
        for block in perm.blocks:
            out_stream.write("{0:+} ".format(block.signed_id()))
        out_stream.write("$\n")
