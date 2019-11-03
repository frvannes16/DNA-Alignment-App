from django.db import models

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
    status = models.CharField(max_length=30, choices=STATUS_CHOICES, default=PROCESSING)
    search_string = models.TextField(null=False, blank=True)
    searches_issued = models.IntegerField(null=False)
    submitted_at = models.DateTimeField(auto_now=True)

class Result(models.Model):
    search = models.ForeignKey(Search, models.CASCADE)
    is_best_match = models.BooleanField(default=False)
    completed_at = models.DateTimeField(auto_now=True)
    match_found = models.BooleanField(default=False)
    match_score = models.DecimalField(null=False, decimal_places=2, max_digits=30)
    match_details = models.TextField(null=False, blank=False)

