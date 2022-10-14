from sys import argv
from xml.etree import ElementTree as et

def edit_PhysiCell_settings_XML_for_global_signal(*args):
  tree = et.parse("./config/PhysiCell_settings.xml")
  keys = args[0][0::2]
  values = args[0][1::2]
  tree.find(".//leaders_have_global_signal").text = 'false'
  tree.find(".//followers_have_global_signal").text = 'false'
  tree.find(".//weight_towards_chemotaxis_cues").text = '0'
  for i in range(len(keys)):
      name = ".//" + keys[i]
      tree.find(name).text = values[i]
      # tree.find('idinfo/timeperd/timeinfo/rngdates/begdate').text = '1/1/2011'
  tree.write("./config/PhysiCell_settings.xml")


if __name__ == "__main__":
  edit_PhysiCell_settings_XML_for_global_signal(argv[1:])
