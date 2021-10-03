
from obspy.core.inventory.inventory import read_inventory

XMLDIR="/cea/data/outils/dase/seiscomp/xml"
network,station,location,channel="XX","AND0","00","DDF"
network,station,location,channel="RD","2HELA","10","LAN"

inv=read_inventory('%s/%s.%s.%s.%s.xml' % (XMLDIR,network,station,location, channel))
resp=inv[0][0][0].response

SENSIB=1.0/inv[0][0][0].response._get_overall_sensitivity_and_gain(frequency=1.0, output="ACC")[1]
print (SENSIB)

print (resp)

resp.plot(0.0000001, output="ACC")



