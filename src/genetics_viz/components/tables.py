"""Table components and slots for genetics-viz."""

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
            <span v-if="col.value === 'present'" style="display: flex; align-items: center; gap: 4px;">
                <q-icon name="check_circle" color="green" size="sm">
                    <q-tooltip>Validated as present</q-tooltip>
                </q-icon>
                <span v-if="props.row.ValidationInheritance === 'de novo' || props.row.Inheritance === 'de novo'" style="font-weight: bold;">dnm</span>
            </span>
            <q-icon v-else-if="col.value === 'absent'" name="cancel" color="red" size="sm">
                <q-tooltip>Validated as absent</q-tooltip>
            </q-icon>
            <q-icon v-else-if="col.value === 'uncertain'" name="help" color="orange" size="sm">
                <q-tooltip>Validation uncertain or different</q-tooltip>
            </q-icon>
            <q-icon v-else-if="col.value === 'conflicting'" name="bolt" color="amber-9" size="sm">
                <q-tooltip>Conflicting validations</q-tooltip>
            </q-icon>
        </template>
        <template v-else>
            {{ col.value }}
        </template>
    </q-td>
</q-tr>
"""
