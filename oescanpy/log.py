import logging
from functools import wraps
from datetime import datetime
import sys
from pprint import pformat

def create_logging_handler(debug: bool) -> logging.Handler:
    """
    Create a handler for logging purposes.

    :param debug: Does debug information need to be logged
    :return: The logging handler.
    """
    # Disable default logging handlers
    logging.basicConfig(handlers=[], level=logging.WARNING)

    # Create our custom handler
    ch = logging.StreamHandler(stream=sys.stderr)
    ch.setLevel(logging.DEBUG if debug else logging.INFO)
    ch.setFormatter(
        logging.Formatter("\n%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    )
    return ch


def log_function_call(func):
    """
    A decorator to measure and log function execution time and arguments.

    Args:
        func (function): The function to be decorated.

    Returns:
        function: The wrapped function with timing and args logging capability.
    """

    @wraps(func)
    def wrapper(*args, **kwargs):
        start_time = datetime.now()

        try:
            return func(*args, **kwargs)
        finally:
            elapsed_time = (datetime.now() - start_time).total_seconds()
            logging.info(f"函数 {func.__name__} 执行耗时: {elapsed_time:.4f}秒")
    return wrapper
