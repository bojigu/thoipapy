
def get_set_lists(s: dict):

    if s["set_list"] == True:
        set_list = ["1"]
    elif s["set_list"] == False:
        set_list = ["0"]
    elif isinstance(s["set_list"], int):
        set_list = [str(s["set_list"])]
    elif isinstance(s["set_list"], str):
        set_list = s["set_list"].split(",")
    else:
        raise ValueError("set_list type is not correct {} ({})".format(s["set_list"], type(s["set_list"])))

    return set_list