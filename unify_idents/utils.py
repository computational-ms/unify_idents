"""Collection of utils."""


def merge_and_join_dicts(list_of_dicts, delimiter):
    """Merge list of dicts with identical keys as strings into single merged dict.

    Args:
        list_of_dicts (list): list of dicts which are to be merged
        delimiter (str): delimiter

    Returns:
        dict: preserved original keys with all values merged as strings with delimiter

    """
    return {
        key: delimiter.join([str(d.get(key)) for d in list_of_dicts])
        for key in set().union(*list_of_dicts)
    }
