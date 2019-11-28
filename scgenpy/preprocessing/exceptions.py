class UnboundAttributeError(UnboundLocalError):
    """
    Raised when a reference is made to a class attribute that is either not initialised or
    initialised as None, This is a subclass of UnboundLocalError.
    """
