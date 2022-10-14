from sys import argv
from xml.etree import ElementTree as et

def edit_PhysiCell_settings_XML(*args):
  tree = et.parse("./config/PhysiCell_settings.xml")
  keys = args[0][0::2]
  values = args[0][1::2]
  tree.find(".//max_time").text = '720'
  tree.find(".//fibronectin_strip_initial_condition").text = 'false'
  tree.find(".//fibronectin_kernel_max_value").text = '10.0'
  for i in range(len(keys)):
      name = ".//" + keys[i]
      tree.find(name).text = values[i]
      if name == ".//followers_cell_cell_repulsion_strength" and values[0] not in ['30', '32', '35', '38', '39', '42', '43', '46', '47', '51', '52', '55', '56','61','62','63','64','65','66','67','68','69','71','72','73','75','76','77']:
          tree.find(".//leaders_cell_cell_repulsion_strength").text = values[i]
      elif name == ".//followers_cell_cell_repulsion_strength" and values[0] in ['30', '35', '39', '43', '47', '52', '56','62','65','68', '72','76']:
         tree.find(".//leaders_cell_cell_repulsion_strength").text = '2'
      elif name == ".//followers_cell_cell_repulsion_strength" and values[0] in ['32', '38', '42', '46', '51', '55','61','64','67','71','75']:
         tree.find(".//leaders_cell_cell_repulsion_strength").text = '0.05'
      elif name == ".//followers_cell_cell_repulsion_strength" and values[0] in ['63','66','69','73','77']:
            tree.find(".//leaders_cell_cell_repulsion_strength").text = '0.5'
      if name == ".//scenario" and values[i] == '12':
          tree.find(".//max_time").text = '1440'
      if name == ".//scenario" and int(values[i]) > 21 and int(values[i]) < 29:
          tree.find(".//fibronectin_strip_initial_condition").text = 'true'
      if name == ".//scenario" and values[i] in ['31', '36', '48', '57', '60', '70', '71', '72', '73']:
          tree.find(".//fibronectin_kernel_max_value").text = '100.0'
      # tree.find('idinfo/timeperd/timeinfo/rngdates/begdate').text = '1/1/2011'
  tree.write("./config/PhysiCell_settings.xml")


if __name__ == "__main__":
  edit_PhysiCell_settings_XML(argv[1:])
