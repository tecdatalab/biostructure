import numpy as np
import torch
from torch.nn import Module
from torch.nn.functional import cross_entropy, softmax


class CustomLoss(Module):
    def __init__(self, name, num_classes, device, alpha, beta, f16=False):
        """
        A wrapper Module for a custom loss function
        """
        super(CustomLoss, self).__init__()
        self.name = name
        self.num_classes = num_classes
        self.alpha = alpha
        self.beta = beta
        self.device = device
        self.half_precision = f16

    def tversky_loss(self, pred, target, alpha, beta, gamma=0.75):
        """
        Calculate the Tversky loss for the input batches
        :param pred: predicted batch from model
        :param target: target batch from input
        :param alpha: multiplier for false positives
        :param beta: multiplier for false negatives
        :return: Tversky loss
        """
        target_oh = torch.eye(self.num_classes)[target.squeeze(1)]
        target_oh = target_oh.permute(0,4,1,2,3).float()
        probs = softmax(pred, dim=1)
        target_oh = target_oh.type(pred.type())
        dims = (0,) + tuple(range(2, target.ndimension()))
        inter = torch.sum(probs * target_oh, dims)
        fps = torch.sum(probs * (1 - target_oh), dims)
        fns = torch.sum((1 - probs) * target_oh, dims)
        t = (inter / (inter + (alpha * fps) + (beta * fns))).mean()
        return (1 - t)**gamma

    def dice_loss(self, pred, target):
        """
        Calculate the Tversky loss for the input batches
        :param pred: predicted batch from model
        :param target: target batch from input
        :param alpha: multiplier for false positives
        :param beta: multiplier for false negatives
        :return: Tversky loss
        """
        target_oh = torch.eye(self.num_classes, device=self.device)[target.squeeze(1)]
        target_oh = target_oh.permute(0,4,1,2,3).float()
        probs = softmax(pred, dim=1)
        target_oh = target_oh.type(pred.type())
        dims = (0,) + tuple(range(2, target.ndimension()))
        inter = torch.sum(probs * target_oh, dims)[1:]
        fps = torch.sum(probs * (1 - target_oh), dims)[1:]
        fns = torch.sum((1 - probs) * target_oh, dims)[1:]
        t = (inter / (inter + (0.5 * fps) + (0.5 * fns))).mean()
        return 1 - t

    def class_dice(self, pred, target, epsilon=1e-6):
        """
        Calculate DICE coefficent for each class
        :param pred: predicted batch from model
        :param target: target batch from input
        :param epsilon: very small number to prevent divide by 0 errors
        :return: list of DICE loss for each class
        """
        pred_class = torch.argmax(pred, dim=1)
        dice = np.ones(self.num_classes)
        for c in range(self.num_classes):
            p = (pred_class == c)
            t = (target == c)
            inter = (p * t).sum().float()
            union = p.sum() + t.sum() + epsilon
            d = 2 * inter / union
            dice[c] = 1 - d
        if self.half_precision:
            return torch.from_numpy(dice).half()
        else:
            return torch.from_numpy(dice).float()

    def unified_loss(self, pred, target, gamma, sigma, lambda_)
        ce = cross_entropy(pred, target, weight=self.class_dice(pred,target).to(self.device), reduction='none')
        tv = self.tversky_loss(pred, target, sigma, 1-sigma, gamma)
        probs = softmax(pred, dim=1)
        dims = (0,) + tuple(range(2, target.ndimension()))
        print(ce.size())
        print(probs.size())
        ce *= (1 - probs) ** gamma  # focal loss factor
        focal_ce = torch.sum(ce, dim=dims).mean()


        loss = (lambda_ * focal_ce) + ((1-lambda_) * tv)


    def forward(self, pred, target, cross_entropy_weight=0.5,
                tversky_weight=0.5):
        """
        Calculate the custom loss
        :param pred: predicted batch from model
        :param target: target batch from input
        :param cross_entropy_weight: multiplier for cross entropy loss
        :param tversky_weight: multiplier for tversky loss
        :return: loss value for batch
        """
        if self.name == 'CrossEntropy':
            loss = cross_entropy(pred, target,
                               weight=self.class_dice(pred,target).to(self.device))
            return loss
        elif self.name == 'Dice':
            loss = self.dice_loss(pred, target)
            return loss
 
        elif self.name == 'Tversky':
            loss = self.tversky_loss(pred, target, self.alpha, self.beta)
            return loss
        elif self.name == 'Unified':
            loss = self.unified_loss(pred, target, self.gamma, self.sigma, self.lambda_)
            return loss
