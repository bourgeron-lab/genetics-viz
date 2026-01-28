"""Table components and slots for genetics-viz."""

from genetics_viz.utils.gene_scoring import get_gene_scorer

# Initialize gene scorer for badge coloring
_gene_scorer = get_gene_scorer()


def create_gene_badge_data(symbol_or_ensg: str) -> dict:
    """Create badge data for a gene symbol or ENSG ID."""
    color = _gene_scorer.get_gene_color(symbol_or_ensg)
    tooltip = _gene_scorer.get_gene_tooltip(symbol_or_ensg)
    return {
        "label": symbol_or_ensg,
        "color": color,
        "tooltip": tooltip,
    }


# Custom slot for validation table with view button and validation icons
VALIDATION_TABLE_SLOT = r"""
<q-tr :props="props">
    <q-td key="actions" :props="props">
        <q-btn 
            flat 
            dense 
            size="sm" 
            icon="visibility" 
            color="blue"
            @click="$parent.$emit('view_variant', props.row)"
        >
            <q-tooltip>View in IGV</q-tooltip>
        </q-btn>
    </q-td>
    <q-td v-for="col in props.cols.filter(c => c.name !== 'actions')" :key="col.name" :props="props">
        <template v-if="col.name === 'Validation'">
            <span v-if="col.value === 'present' || col.value === 'in phase MNV'" style="display: flex; align-items: center; gap: 4px;">
                <q-icon name="check_circle" color="green" size="sm">
                    <q-tooltip>Validated as {{ col.value }}</q-tooltip>
                </q-icon>
                <span v-if="props.row.ValidationInheritance === 'de novo' || props.row.Inheritance === 'de novo'" style="font-weight: bold;">dnm</span>
                <span v-else-if="props.row.ValidationInheritance === 'homozygous' || props.row.Inheritance === 'homozygous'" style="font-weight: bold;">hom</span>
                <span v-if="col.value === 'in phase MNV'" style="font-size: 0.75em; color: #666;">MNV</span>
            </span>
            <q-icon v-else-if="col.value === 'absent'" name="cancel" color="red" size="sm">
                <q-tooltip>Validated as absent</q-tooltip>
            </q-icon>
            <q-icon v-else-if="col.value === 'uncertain' || col.value === 'different'" name="help" color="orange" size="sm">
                <q-tooltip>Validation uncertain or different</q-tooltip>
            </q-icon>
            <q-icon v-else-if="col.value === 'conflicting'" name="bolt" color="amber-9" size="sm">
                <q-tooltip>Conflicting validations</q-tooltip>
            </q-icon>
        </template>
        <template v-else-if="col.name.toLowerCase().includes('symbol') || col.name.toLowerCase().includes('gene')">
            <template v-if="props.row[col.name + '_badges'] && props.row[col.name + '_badges'].length > 0">
                <q-badge 
                    v-for="(badge, idx) in props.row[col.name + '_badges']" 
                    :key="idx"
                    :label="badge.label" 
                    :style="'background-color: ' + badge.color + '; color: ' + (badge.color === '#ffffff' ? 'black' : 'white') + '; font-size: 0.875em; padding: 4px 8px;'"
                    class="q-mr-xs q-mb-xs"
                >
                    <q-tooltip>{{ badge.tooltip }}</q-tooltip>
                </q-badge>
            </template>
            <template v-else>
                <span>{{ col.value || '-' }}</span>
            </template>
        </template>
        <template v-else>
            {{ col.value }}
        </template>
    </q-td>
</q-tr>
"""
