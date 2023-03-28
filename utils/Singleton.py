#!/usr/bin/env python3
# -*- coding: utf-8 -*-


class Singleton(object):
    _Instance = {}

    def __init__(self, cls):
        self.cls = cls

    def __call__(self, *args, **kwargs):
        # mlogger.debug("call Singleton")
        instance = self._Instance.get(self.cls, None)
        if not instance:
            # mlogger.debug("already exist Singleton")
            instance = self.cls(*args, **kwargs)
            self._Instance[self.cls] = instance
        return instance

    def __getattr__(self, item):
        return getattr(self.cls, item, None)
