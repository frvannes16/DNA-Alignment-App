from __future__ import absolute_import
import os
from celery import Celery
from django.conf import settings

# Set the default Django settings module for celery workers
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'dna_alignment_manager.settings')
app = Celery('dna_alignment_manager')

app.config_from_object('django.conf:settings')
app.autodiscover_tasks(lambda: settings.INSTALLED_APPS)

@app.task(bind=True)
def debug_task(self):
    print('Request: {0!r}'.format(self.request))