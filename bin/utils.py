import logging
import os


def create_logger(logger_name, log_file=False):
    '''
    Creates a custom logger element
    '''
    # format for logging message
    FORMATTER = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'

    # create the logger and set the level to Debug
    logger = logging.getLogger(logger_name)
    logger.setLevel(level=logging.DEBUG)
    # create the handler and set the level to Info
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    # set format
    console_handler.setFormatter(logging.Formatter(FORMATTER))
    # add the handler to the logger
    logger.addHandler(console_handler)

    if log_file:
        file_handler = logging.FileHandler(filename='Logs.log')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(logging.Formatter(FORMATTER))
        logger.addHandler(file_handler)

    return(logger)


def get_basename(file_name):
    '''
    Retrieves the basename of the file to allow flexible output names
    '''
    # retrives the basename of the input file
    if "/" in file_name:
        basename = file_name.split("/")[-1].split(".")[0]
    else:
        basename = file_name.split(".")[0]
    
    return(basename)


def check_output_dir(output_dir):
    """
    Checks if output dir exists, if not create output dir
    Checks that output dir is formatted correctly with '/' on the end - if not, adds '/'
    """
    # checks if output directory exists or creates one
    if not (os.path.exists(output_dir)):
        os.mkdir(output_dir)
        logging.info('Creating output directory')

    if output_dir.endswith("/"):
        output_dir = output_dir
    else:
        output_dir = output_dir + "/"
    # checks if output directory exists after attempting to create one
    if not (os.path.exists(output_dir)):
        logging.info('Output dir could not be created - check permissions & create output dir manually if necessary')
    
    return(output_dir)