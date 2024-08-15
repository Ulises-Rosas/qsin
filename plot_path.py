import numpy as np
import matplotlib.pyplot as plt

errors_path = np.loadtxt("./test_data/linear_batches_testErrors.csv", delimiter=',')
elnet_path  = np.loadtxt("./test_data/linear_batches_elnetPath.csv", delimiter=',')

# errors_path = np.loadtxt("./test_data/non_linear_batches_testErrors.csv", delimiter=',')
# elnet_path  = np.loadtxt("./test_data/non_linear_batches_elnetPath.csv", delimiter=',')
fs = 18
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7))

ax1.plot(elnet_path[:,0], elnet_path[:,1:], marker='o', alpha=0.8)
ax1.set_xscale('log')
ax1.set_ylabel('Coefficient', fontsize=fs)
ax1.axhline(0, color='black', lw=2)
ax1.set_title('Elastic Net Path', fontsize=fs)

ax2.plot(errors_path[:,0], errors_path[:,1], marker='o', alpha=0.8)
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlabel('$\\lambda$', fontsize=fs)
ax2.set_ylabel('RMSE', fontsize=fs)
ax2.set_title('Test Error', fontsize=fs)
plt.tight_layout()
# save the figure
plt.savefig('./imgs/linear_batches_elnetPath.png', dpi=300)

