###############################################################################

def _parse_control(
    control: Mapping[str, Union[float, int]]
) -> Mapping[str, Union[float, int]]:
    if control["precision"] == "float32":
        dtype = torch.float32
    else:
        dtype = torch.float64
    return {
        "max_iter": int(control["max_itr"]),
        "history_size": int(control["history_size"]),
        "tolerance_grad": float(control["tolerance_grad"]),
        "tolerance_change": float(control["tolerance_change"]),
        "dtype": dtype
    }

###############################################################################
