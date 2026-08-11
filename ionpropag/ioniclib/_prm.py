def parse_prs(prs):

    # this must be local, otherwise they are not freed
    class PRMstructure(object):
        name = None
        members = []

    class PRMmember(object):
        pass

    prs = prs.replace("{"," { ").replace("}"," } ").replace("="," = ").replace(","," , ")
    root = None
    curr = None
    parent = None
    next_is_name  = False
    next_is_value = ""
    next_values = None
    for token in prs.split():
        if token == "=": continue
        if token == ",": continue
        if next_is_name:
            curr.name = token
            next_is_name = False
        elif next_is_value != "":
            if token == '{':
                # vector-valued
                next_values = []
            elif token == '}':
                setattr(curr,next_is_value,next_values)
                next_is_value = ""
                next_values = None
            elif next_values is not None:
                next_values.append(token)
            else:
                setattr(curr,next_is_value,token)
                next_is_value = ""
        elif token == "$structure":
            root = PRMstructure()
            curr = root
            next_is_name = True
        elif token == "$member":
            s = PRMmember()
            curr.members.append(s)
            parent = curr
            curr = s
            next_is_name = True
        elif token in ["$type","$default"]:
            next_is_value = token[1:]
        elif token == "}":
            curr = parent

    return root
