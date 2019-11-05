from django.db import models
from django.core.serializers.json import DjangoJSONEncoder
from typing import Dict
import json


class Search(models.Model):
    PROCESSING = 'PR'
    MATCH_FOUND = 'MF'
    NO_MATCH = 'NM'
    FAILED = 'FL'
    STATUS_CHOICES = (
        (PROCESSING, 'PROCESSING'),
        (MATCH_FOUND, 'MATCH_FOUND'),
        (NO_MATCH, 'NO_MATCH'),
        (FAILED, 'FAILED')
    )

    status_dict = {choice: value for (choice, value) in STATUS_CHOICES}

    status = models.CharField(
        max_length=30, choices=STATUS_CHOICES, default=PROCESSING)
    search_string = models.TextField(null=False, blank=True)
    searches_issued = models.IntegerField(null=False)
    submitted_at = models.DateTimeField(auto_now=True)

    def to_json_dict(self) -> Dict:
        data = {
            "searchId": self.id,
            "status": self.status_dict[self.status],
            "searchString": self.search_string
        }
        best_match = self.result_set.filter(
            match_found=True, is_best_match=True).first()
        if best_match is not None:
            data['match'] = best_match.to_json_dict()
        return data


class Result(models.Model):
    search = models.ForeignKey(Search, models.CASCADE)
    is_best_match = models.BooleanField(default=False)
    completed_at = models.DateTimeField(auto_now=True)
    match_found = models.BooleanField(default=False)
    match_score = models.DecimalField(
        null=False, decimal_places=2, max_digits=30)
    protein = models.CharField(max_length=500, null=True, blank=False)
    match_start = models.IntegerField(null=True)
    match_end = models.IntegerField(null=True)

    def to_json_dict(self) -> Dict:
        return {
            'protein': self.protein,
            'startPos': self.match_start,
            'endPos': self.match_end
        }
