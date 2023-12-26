README for the subjective performance evaluation code

These code are used to implement the subjective performance evaluation.
Authors: Gongping Huang
Date: 06/16/2017

USAGE: source = se_evaluate_score(ref_signal, target_signal);

Inputs:
ref_signal:       the reference signal (can be filename)
target_signal:    the target signal (can be filename)
do_sync:          choose do synchronization or not
----------------------------------------------------------------
Output:
source:           the subjective performance score
----------------------------------------------------------------
Notice:

1) The reference signal and target_signal in the input can be
either a real value vector or a filename.

2) In implementing the evaluation, the ref_signal and the
target_signal should be synchronized/aligned, however,this
function can help synchronize the two signals.

3) The synchronization is achieved with a time delay estimation
with the well known GCC-PHAT method, however, the accuracy of the
TDE cannot be guaranteed. Consequently, it is better to synchronize
the reference signal and target_signal before call this function.
In this case, you can just set do_sync = 0.


Reference:
[1] Kinoshita, Keisuke, et al. "A summary of the REVERB challenge:
state-of-the-art and remaining challenges in reverberant speech
processing research." EURASIP Journal on Advances in Signal
Processing 2016.1 (2016): 1-19.

Authors: Gongping Huang
Date: 06/16/2017