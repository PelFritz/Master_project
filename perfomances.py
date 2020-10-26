"""This is script can be used to calculate the precision, sensitivty, specificity and f1 score"""


def sensitivity(tp, fn):
    return tp/(tp + fn)


def specificty(tn, fp):
    return tn/(tn + fp)


def precision(tp, fp):
    return tp/(tp + fp)


def f1_score(tp, fp, fn):
    prec = precision(tp, fp)
    recall = sensitivity(tp, fn)
    return 2*(prec*recall/(prec + recall))


tp = 496262
tn = 16122
fp = 88944
fn = 40096

print('sensitivty:\t', sensitivity(tp, fn))
print('Specifity:\t', specificty(tn, fp))
print('Precision:\t', precision(tp, fp))
print('F1 score:\t', f1_score(tp, fp, fn))
