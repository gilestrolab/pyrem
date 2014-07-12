import numpy as np
import cpython_proto as cp
random_walk = np.cumsum(np.random.normal(0,1,(int(2e3))))
print "?"
for i in range(20):
    print cp.ap_entropy(random_walk, 2,1.5)
