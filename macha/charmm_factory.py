import natsort


class CharmmFactory:
    def __init__():
        pass

    def createHeader(segids, env):

        header = ""
        for segi in natsort.natsorted(segids):
            if segi.startswith("RNA"):
                header += f"""
! Read {segi.upper()}
open read card unit 10 name {segi.lower()}.crd
read sequence coor card unit 10 resid
generate {segi.upper()} setup warn first 5TER last 3TER

open read unit 10 card name {segi.lower()}.crd
read coor unit 10 card resid
        """
            elif env != "waterbox" and segi not in [
                "HETA",
                "IONS",
                "SOLV",
                "WATA",
                "WATB",
                "WATC",
            ]:
                header += f"""
! Read {segi.upper()}
open read card unit 10 name {segi.lower()}.crd
read sequence coor card unit 10 resid
generate {segi.upper()} setup warn first 5TER last 3TER

open read unit 10 card name {segi.lower()}.crd
read coor unit 10 card resid
        """
            elif segi == "HETA":
                header += f"""
bomlev -1  ! not ideal but necessary for taking three membered rings into account
! Read {segi.upper()}
open read card unit 10 name {segi.lower()}.crd
read sequence coor card unit 10 resid
generate {segi.lower()} setup warn first none last none

open read unit 10 card name {segi.lower()}.crd
read coor unit 10 card resid
bomlev 0
        """
            else:
                pass

        return header
