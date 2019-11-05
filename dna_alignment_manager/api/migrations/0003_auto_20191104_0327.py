# Generated by Django 2.2.6 on 2019-11-04 03:27

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('api', '0002_result_match_details'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='result',
            name='match_details',
        ),
        migrations.AddField(
            model_name='result',
            name='match_end',
            field=models.IntegerField(null=True),
        ),
        migrations.AddField(
            model_name='result',
            name='match_start',
            field=models.IntegerField(null=True),
        ),
        migrations.AddField(
            model_name='result',
            name='protein',
            field=models.CharField(max_length=500, null=True),
        ),
    ]
