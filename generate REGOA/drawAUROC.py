import pickle
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot,savefig
from sklearn import svm, datasets
from sklearn.metrics import roc_curve, auc  ###计算roc和auc


def getauroc(y,prob):
    fpr, tpr, threshold = roc_curve(y, prob)  ###计算真正率和假正率
    roc_auc = auc(fpr, tpr)  ###计算auc的值
    lw = 1
    plt.plot(fpr, tpr, color='darkorange',
             lw=lw, label='ROC curve (area = %0.3f)' % roc_auc)  ###假正率为横坐标，真正率为纵坐标做曲线
    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC of metapath2vec in GO&GRN_'+Function)
    plt.legend(loc="lower right")
    return roc_auc

if __name__=='__main__':
    Function = 'MF'
    plt.figure(figsize=(10, 10))
    g = open('./datahg/AUROCs/' + Function + '_standard_score_GOGRNPPI.pkl', 'rb+')
    scores = pickle.load(g)
    fpr, tpr, threshold = roc_curve(scores[0], scores[1])
    roc_auc = auc(fpr, tpr)
    lw = 2
    plt.plot(fpr, tpr, color='#FF0000',
             lw=lw, label='RE-GOA (GRN+PPI+GO) \n(AUROC = %0.3f)' % roc_auc)  ###假正率为横坐标，真正率为纵坐标做曲线

    g = open('./datahg/AUROCs/' + Function + '_standard_score_GOGRN.pkl', 'rb+')
    scores = pickle.load(g)
    fpr, tpr, threshold = roc_curve(scores[0], scores[1])
    roc_auc = auc(fpr, tpr)

    plt.plot(fpr, tpr, color='#00FF00',
             lw=lw, label='GO+GRN (AUROC = %0.3f)' % roc_auc)  ###假正率为横坐标，真正率为纵坐标做曲线

    g = open('./datahg/AUROCs/' + Function + '_standard_score_GOPPI.pkl', 'rb+')
    scores = pickle.load(g)
    fpr, tpr, threshold = roc_curve(scores[0], scores[1])
    roc_auc = auc(fpr, tpr)

    plt.plot(fpr, tpr, color='#0000FF',
             lw=lw, label='GO+PPI (AUROC = %0.3f)' % roc_auc)  ###假正率为横坐标，真正率为纵坐标做曲线

    g = open('./datahg/diffpaths/' + Function + '_standard_score-noRT.pkl', 'rb+')
    scores = pickle.load(g)
    fpr, tpr, threshold = roc_curve(scores[0], scores[1])
    roc_auc = auc(fpr, tpr)

    plt.plot(fpr, tpr, color='orange',
             lw=lw, label='RE-GOA without RE-terms \n(AUROC = %0.3f)' % roc_auc)  ###假正率为横坐标，真正率为纵坐标做曲线

    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.tick_params(labelsize=25)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.0])
    plt.xlabel('False Positive Rate',fontsize=30)
    plt.ylabel('True Positive Rate',fontsize=30)
    plt.title('ROC of Gene Function Prediction in ' + Function+'\n',fontsize=30)
    plt.legend(loc="lower right",fontsize=25)
    #plt.show()
    plt.savefig("./datahg/AUROCs/AUROCs_" + Function + ".pdf")

