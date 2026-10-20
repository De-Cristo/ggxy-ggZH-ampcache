"""Decay-channel descriptions for exported models; never change weights or map data."""
from copy import deepcopy


def channel_card(card, channel):
    if channel not in ("ll", "vv"):
        raise ValueError(channel)
    result = deepcopy(card)
    if channel == "vv":
        result["channel"] = "Z -> nu nubar, summed over three neutrino flavours"
        result["source"] = result["source"].replace("l+l-H", "nu nubar H").replace("m_ll", "m_nunu")
        suffix = "; Z -> nu nubar, three flavours (set vv13600)"
        if not result["target"].endswith(suffix):
            result["target"] += suffix
        result["objects"] = result["objects"].replace("decay leptons' parent", "decay neutrinos' parent")
        result["classifier"]["weight"] = result["classifier"]["weight"].replace("m_ll", "m_nunu")
    return result


def channel_map(cset, channel):
    if channel not in ("ll", "vv"):
        raise ValueError(channel)
    result = deepcopy(cset)
    if channel == "vv":
        for correction in result["corrections"]:
            correction["description"] = correction["description"].replace("Z(ll)", "Z(nu nubar)").replace("m_ll", "m_nunu")
    return result
