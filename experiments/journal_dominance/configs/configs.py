#!/usr/bin/env python
# encoding: utf-8

class Config:    
    def __init__(self, folder, nick, config, type):
        self.folder = folder
        self.nick = nick
        self.config = ["--search", config]
        self.type = type

    def __repr__(self):
        return ", ".join(map(str, [self.folder, self.nick, self.config, self.type]))



import journal_1_analysis_dominance_relations

CONFIGS = {"journal1" : journal_1_analysis_dominance_relations.CONFIGS}


for conf in CONFIGS:
    print conf
    print ""
