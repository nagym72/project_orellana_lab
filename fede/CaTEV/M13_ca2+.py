cmd.remove ("resn CLA or resn POT or resn TIP3")
cmd.show ("spheres", "resn CAL")
cmd.color ("yellow", "resn CAL")

cmd.select("CaM", "resi 65:231")
cmd.color("red", "CaM")

cmd.select("protease", "resi :62 or resi 235:402")
cmd.color("skyblue", "protease")

cmd.select("MK2", "resi 144:160")
cmd.color("yelloworange", "MK2")

cmd.select("linkers", "resi 63+64+142+143+161+162+232+233+234")
cmd.color("forest", "linkers")

cmd.select("EF1", "resi 82+84+86+93")
cmd.color("green", "EF1")

cmd.select("EF2", "resi 118+120+122+129")
cmd.color("cyan", "EF2")

cmd.select("EF3", "resi 176+178+187")
cmd.color("cyan", "EF3")

cmd.select("EF4", "resi 212+214+216+223")
cmd.color("green", "EF4")

cmd.deselect()