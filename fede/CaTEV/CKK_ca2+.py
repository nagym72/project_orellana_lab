cmd.remove ("resn CLA or resn POT or resn TIP3")
cmd.show ("spheres", "resn CAL")
cmd.color ("yellow", "resn CAL")

cmd.select("CaM", "resi 64:240")
cmd.color("red", "CaM")

cmd.select("protease", "resi :61 or resi 244:401")
cmd.color("skyblue", "protease")

cmd.select("MK2", "resi 144:169")
cmd.color("yelloworange", "MK2")

cmd.select("linkers", "resi 62+63+142+143+169+170+171+241+242+243")
cmd.color("forest", "linkers")

cmd.select("EF1", "resi 82+84+86+93")
cmd.color("green", "EF1")

cmd.select("EF2", "resi 118+120+122+129")
cmd.color("cyan", "EF2")

cmd.select("EF3", "resi 185+187+196")
cmd.color("cyan", "EF3")

cmd.select("EF4", "resi 221+223+225+232")
cmd.color("green", "EF4")

cmd.deselect()