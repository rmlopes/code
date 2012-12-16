from logging.handlers import RotatingFileHandler

class WriteRotateFileHandler(RotatingFileHandler):
        def __init__(self, filename):
                RotatingFileHandler.__init__(self, filename, 'w',0,0)

        def emit(self, record):
            self.doRollover()
            RotatingFileHandler.emit(self, record)


#if __name__ == '__main__':
       # wrfh = WriteRotateFileHandler('../../results/code_circuit_0.save')
