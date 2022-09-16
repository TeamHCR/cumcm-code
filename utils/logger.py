import logging
from colorlog import ColoredFormatter


def get_logger(name=__name__):
    logger_base = logging.getLogger(name)
    logger_base.setLevel(logging.DEBUG)
    stream_handler = logging.StreamHandler()
    color_formatter = ColoredFormatter(
        # '%(log_color)s[%(pathname)s:%(lineno)s %(levelname)s] %(message)s'
        '%(log_color)s[%(funcName)s %(levelname)s] %(message)s'
        )
    stream_handler.setFormatter(color_formatter)
    logger_base.addHandler(stream_handler)
    return logger_base


logger = get_logger(__name__)

if __name__ == '__main__':
    logger.debug('debug message')
    logger.info('info message')
    logger.warn('warn message')
    logger.error('error message')
    logger.critical('critical message')
