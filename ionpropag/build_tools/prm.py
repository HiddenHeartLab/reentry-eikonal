import json
import re

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

def generate_parameters(prs_file, json_file):
    print(f"Parsing {prs_file} parameters...")

    prs = parse_prs(prs_file.read_text())

    params = {
        "name": prs.name,
        "members": [],
    }

    for member in prs.members:
        val = member.default
        typ = member.type

        match = re.match(
            r"^\$(float|double|int|boolean)(?:|\[([0-9]+)\])$",
            typ,
        )

        if not match:
            raise RuntimeError(
                f"Cannot parse {member.type} in PRM file"
            )

        typ = match.group(1)

        if typ == "int":
            conv = int
        elif typ in ("float", "double"):
            conv = float
        elif typ == "boolean":
            conv = lambda x: 1 if x == "TRUE" else 0

        if isinstance(val, list):
            val = [conv(v) for v in val]
        else:
            val = conv(val)

        arr = int(match.group(2)) if match.group(2) else 0

        params["members"].append(
            (member.name, typ, arr, val)
        )

    # PRM reverses parameter order in the generated .h file.
    params["members"].reverse()

    json_file.write_text(
        json.dumps(params, indent=2) + "\n"
    )
