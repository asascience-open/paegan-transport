import logging

try:
    PROGRESS=15
    logging.PROGRESS = PROGRESS
    logging.addLevelName(PROGRESS, 'PROGRESS')
    def progress(self, message, *args, **kws):
        if self.isEnabledFor(PROGRESS):
            self._log(PROGRESS, message, args, **kws)
    logging.Logger.progress = progress
except:
    pass

__version__ = '0.2'